function [Bo,Fo,To,Fco,Tco] = raspecgram(varargin)
%   raspecgram - compute a reassigned spectrogram
%   of a signal using the method of Auger and 
%   Flandrin to compute the partial derivatives 
%   from STFTs.
% ---------------------------------------------------
%   [B,F,T] = raspecgram(A,NFFT,Fs,window,olap)
%
%   A =      input signal (samples)
%   NFFT =   FFT size (default min of 256 and length of x)
%   Fs =     sample rate (default 2)
%   window = analysis window (samples, or hannning window 
%            length, default nfft)
%   olap =   overlap of analysis windows (number of samples,
%            default half the window length)
%
%   B =      STFT matrix, columns are discrete Fourier 
%            transforms of short-time slices of the signal.
%   F =      frequency matrix containing reassigned 
%            frequencies (in Hz) for each point in B
%   T =      time matrix containing reassigned times
%            (in seconds) for each point in B
%
%   [B,Fnom,Tnom,Fcor,Tcor] = raspecgram(...)
%
%   Fnom =   matrix containing nominal (not reassigned)
%            frequencies (in Hz) for each element in B
%            (values in a row are all identical)
%   Tnom =   matrix containing nominal (not reassigned)
%            times (in seconds) for each element in B
%            (values in a column are all identical)
%   Fcor =   matrix containing frequency corrections 
%            (in Hz) for each point in B
%   Tcor =   matrix containing time corrections
%            (in seconds) for each point in B
%
%   raspecgram(...)
%
%   If no return arguments are specified, a log magnitude
%   image plot is constructed and displayed.
%
%   See also crossspecgram, findifspecgram, raspectrum.

%---------------------------------------------------------------------
%	Copyright (c) 2005,2006,2007 by Kelly Fitz and Sean Fulop.
%
%	Please send comments or bug reports to the authors at
%
%       kfitz@cerlsoundgroup.org
%       sfulop@csufresno.edu
%
%   or visit http://www.cerlsoundgroup.org/Kelly/timefrequency.html
%
%   Last updated 20 March 2007 by Fitz.
%
%   This program is free software; you can redistribute it and/or
%   modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation; either version 2
%   of the License, or any later version.
%
%   If you cannot abide the terms of the GNU General Public License,
%   contact the authors regarding an alternate license.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%---------------------------------------------------------------------

% -- check arguments --
% use the auxiliary function, below
error( nargchk(1,5,nargin) )
[msg,A,NFFT,Fs,W,olap] = check_args( varargin );
error( msg )

% --------------- interesting stuff begins here ---------------

% -- compute spectra --
X = specgram(A,NFFT,Fs,W,olap);
magsqrd = abs(X).^2; % used to compute time and freq corrections
nonzero = find(magsqrd>0);
Nw = length(W);

% -- construct reassignment windows --
Wdt = timederivwindow( W, Fs );
Wt = timerampwindow( W, Fs );

% -- compute auxiliary spectra --
Xdt = specgram(A,NFFT,Fs,Wdt,olap);
Xt = specgram(A,NFFT,Fs,Wt,olap);

% -- compute frequency corrections --
[rows,cols] = size(X);
fcorrect = zeros(rows,cols);
fcorrect(nonzero) = -imag(Xdt(nonzero).*conj(X(nonzero)))./magsqrd(nonzero); 

Fk = ((0:rows-1)'*Fs/NFFT)*ones(1,cols);
        % analysis bin frequencies in Hz
        
F = Fk + fcorrect;  % in Hz

% -- compute time corrections --
tcorrect = zeros(rows,cols);
tcorrect(nonzero) = real(Xt(nonzero).*conj(X(nonzero)))./ magsqrd(nonzero);

frametimes = ( ((Nw-1)/2) + (ones(rows,1)*(0:cols-1))*(Nw-olap) ) / Fs;
        % analysis frame times in seconds
        % half the window length ((Nw-1)/2) plus
		% the number (0:cols-1) of frame advances 
		% (Nw-olap) corresponding to the frame, 
		% times the sample period (1/Fs).
T = frametimes + tcorrect; % in seconds

% --------------- interesting stuff ends here ---------------


% -- assign return values, or make a plot --
%
if nargout == 0
    
    % using ratoimage, after cleaning up the data
    nhoriz = max( 500, size(X,2)*2 );
    nvert = max( 400, size(X,1)*2 );
    
    Fmax = .3*Fs; 
    Fmax = 4e6;
    Fmin = 0;
    Tmax = length(A) / Fs;
    Tmin = 0;
    ampthresh = -80;
    
    [X,F,T] = raframe( X, F, T, ampthresh, Fmin, Fmax, Tmin, Tmax );        
    ratoimage( X, F, T, nhoriz, nvert );
    
elseif nargout < 4
    Bo = X;
    Fo = F;
    To = T;
else
    Bo = X;
    Fo = Fk;
    To = frametimes;
    Fco = fcorrect;
    Tco = tcorrect;
end

% --------------- argument checking function ---------------

function [msg,x,nfft,Fs,win,nolap] = check_args(args)
%   Extract and check up to five arguments to raspecgram
%   from an array.
%
%   args =  array of arguments
%   x =     input signal (samples)
%   nfft =  FFT size (default min of 256 and length of x)
%   Fs =    sample rate (default 2)
%   win =   analysis window (samples, or hannning window 
%           length, default nfft)
%   nolap = overlap of analysis windows (number of samples,
%           default half the window length)
%
msg = [];
 
x = args{1};
if (length(args) > 1) && ~isempty(args{2})
    nfft = args{2};
else
    nfft = min(length(x),256);
end
if (length(args) > 2) && ~isempty(args{3})
    Fs = args{3};
else
    Fs = 2;
end
if length(args) > 3 && ~isempty(args{4})
    win = args{4};
else
    if length(nfft) == 1
        win = hanning(nfft);
    else
        msg = 'raspecgram: You must specify a window function.';
    end
end
if length(win) == 1
    win = hanning(win); 
end
if (length(args) > 4) && ~isempty(args{5})
    nolap = args{5};
else
    nolap = ceil(length(win)/2);
end
 
% do error checking
if (length(nfft)==1) && (nfft<length(win)),
    msg = 'raspecgram: Window length cannot be greater than the FFT length.';
end
if (nolap >= length(win)),
    msg = 'raspecgram: Overlap must be less than the window length.';
end
if (length(nfft)==1) && (nfft ~= abs(round(nfft)))
    msg = 'raspecgram: FFT length must be a positive integer.';
end
if (nolap ~= abs(round(nolap))),
    msg = 'raspecgram: Overlap must be a positive integer.';
end
if min(size(x))~=1,
    msg = 'raspecgram: Input signal must be a vector (either row or column).';
end
