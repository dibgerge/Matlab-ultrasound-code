function [mpd, Xo, Fo, To] = mixedpd(varargin)
%   mixedpd - compute the mixed partial derivative of phase
% ---------------------------------------------------------
%   mpd = mixedpd(A,NFFT,Fs,window,olap)
%   [mpd,X] = mixedpd(A,NFFT,Fs,window,olap)
%   [mpd,X,F] = mixedpd(A,NFFT,Fs,window,olap)
%   [mpd,X,F,T] = mixedpd(A,NFFT,Fs,window,olap)
%
%   A = signal
%   NFFT = number of points in spectrum computation
%   Fs = sample rate
%   window = window function for spectrum computation
%   olap =   overlap of analysis windows (number of samples,
%            default half the window length)
%
%   mpd = mixed partial derivatives of phase
%
%   X = STFT matrix, columns are discrete Fourier 
%      	transforms of short-time slices of the signal.
%   F = frequency matrix containing reassigned 
%       frequencies (in Hz) for each point in X
%   T = time matrix containing reassigned times
%       (in seconds) for each point in X
%
%   This gives approximately mpd = -1 in the vicinity
%   of spectral peaks due to sinusoids, and something
%   near mpd = 0 for peaks due to impulses.
%

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
%   Last updated 18 Apr 2007 by Fitz.
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

% -- compute spectra --
X = specgram(A,NFFT,Fs,W,olap);
magsqrd = abs(X).^2; % used to compute time and freq corrections
Nw = length(W);

% -- construct reassignment windows --
WD = timederivwindow(W,Fs);
WT = timerampwindow(W,Fs);
WTD = timerampwindow( timederivwindow(W,Fs), Fs );
    % can use the same ramps as the reassigned spectrogram,
    % need only to scale results by 2 pi (see below)

% -- compute auxiliary spectra --
XD = specgram(A,NFFT,Fs,WD,olap);
XT = specgram(A,NFFT,Fs,WT,olap);
XTD = specgram(A,NFFT,Fs,WTD,olap);

% -- compute mixed partial derivatives --

% add eps to avoid divide by zero?
Xeps = max( X, eps );
term1 = real(XTD./Xeps);
term2 = - real(XD./Xeps) .* real(XT./Xeps);
term3 = imag(XD./Xeps) .* imag(XT./Xeps);
mpd = 2 * pi * (term1 + term2 + term3);
    % scale by 2 pi to normalize (so that sinusoids
    % give values around -1)

if 1 < nargout
    Xo = X;
end

if 2 < nargout
    % -- compute frequency corrections --
    [rows,cols] = size(X);
    magsqrd = abs(Xeps).^2;
    fcorrect = -imag(XD.*conj(Xeps))./magsqrd; 

    Fk = ((0:NFFT/2)'*Fs/NFFT)*ones(1,cols);
            % analysis bin frequencies in Hz

    Fo = Fk + fcorrect;  % in Hz
end

if 3 < nargout
    % -- compute time corrections --
    tcorrect = real(XT.*conj(Xeps))./ magsqrd;

    frametimes = ( ((Nw-1)/2) + (ones(rows,1)*(0:cols-1))*(Nw-olap) ) / Fs;
            % analysis frame times in seconds
            % half the window length ((Nw-1)/2) plus
            % the number (0:cols-1) of frame advances 
            % (Nw-olap) corresponding to the frame, 
            % times the sample period (1/Fs).
    To = frametimes + tcorrect; % in seconds   
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
        msg = 'mixedpd: You must specify a window function.';
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
    msg = 'mixedpd: Window length cannot be greater than the FFT length.';
end
if (nolap >= length(win)),
    msg = 'mixedpd: Overlap must be less than the window length.';
end
if (length(nfft)==1) && (nfft ~= abs(round(nfft)))
    msg = 'mixedpd: FFT length must be a positive integer.';
end
if (nolap ~= abs(round(nolap))),
    msg = 'mixedpd: Overlap must be a positive integer.';
end
if min(size(x))~=1,
    msg = 'mixedpd: Input signal must be a vector (either row or column).';
end
