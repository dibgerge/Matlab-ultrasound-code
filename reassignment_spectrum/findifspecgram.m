function [Bo,Fo,To] = findifspecgram(varargin)
%   findifspecgram - compute a reassigned spectrogram
%   of a signal using finite differences of phase 
%   to approximate the partial derivatives.
% -------------------------------------------------------------------
%   [B,F,T] = findifspecgram(A,NFFT,Fs,window,olap)
%   [B,F] = findifspecgram(A,NFFT,Fs,window,olap)
%   B = findifspecgram(A,NFFT,Fs,window,olap)
%   findifspecgram(A,NFFT,Fs,window,olap)
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
%   See also findifspectrum, raspecgram, raspectrum.

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
%   Last updated 22 Feb 2007 by Fitz.
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
% use the auxiliary function, below, 
error( nargchk(1,5,nargin) )
[msg,A,NFFT,Fs,W,olap] = check_args( varargin );
error( msg )

% pad A with zeros if necessary 
if length(A) <= length(W)
    A = [ A ; zeros( 1 + length(W) - length(A), 1 ) ];
end

% specgram is not convenient for computing cross-spectra,
% manually build up the short-time segment matrices and 
% pass them to fft.
hop = length(W)-olap;
offset = [0 : hop : length(A)-length(W)-1 ];
WW = W*ones(1,length(offset));

% create a delayed-by-one version of the input sample vector
nsamps_win = length(W);
idx = ((1:nsamps_win)'*ones(1,length(offset))) + (ones(nsamps_win,1)*offset);
S = A(idx+1).*WW;
Sdel = A(idx).*WW;

% compute short-time transform surfaces
% (the nominal frame times are used below to compute 
% the corrected times)
STFT = fft(S, NFFT);
STFTdel = fft(Sdel, NFFT);

% construct the frequency-delayed-by-one spectrum
% (delay by rotating, since the spectrum is periodic)
STFTfreqdel = [ STFT(NFFT,:); STFT(1:NFFT-1,:) ];

% -- use only samples below Fs/2 --
Fnominal = ((0:NFFT-1)*Fs/NFFT);
f_select = find(Fnominal < Fs/2);
STFT = STFT(f_select,:);
STFTdel = STFTdel(f_select,:);
STFTfreqdel = STFTfreqdel(f_select,:);
 
% -- approximate derivatives with finite differences --
% compute channelized instaneous frequency, a vector of 
% remapped frequencies
%
% The mod 2pi wraps negative differences up to positive
% differences (the absolute value of the difference can
% never exceed 2pi, since angle has a range of -pi,pi).
argC = mod( angle( STFT ) - angle( STFTdel ), 2*pi );
CIF = (Fs/(2.*pi)) * argC;

% compute local group delay, a vector of time adjustments,
% one for each discrete frequency step determined 
% by fft length
% normalize the local group delay by the
% FFT frequency bin size in radians.
%
% The mod 2pi wraps negative differences up to positive
% differences (the absolute value of the difference can
% never exceed 2pi, since angle has a range of -pi,pi).
%
% Use mod -2pi, because these differences are expected to
% be negative, on the range -2pi,0.
argL = mod( angle( STFT ) - angle( STFTfreqdel ), -2*pi );
LGD = -(NFFT/(Fs*2*pi)) * argL;

% compute the corrected times
[nrows,ncols] = size(STFT);
frames = offset + 1;
        % analysis frame times in samples
Tnominal = ones(nrows,1) * frames/Fs;

X = STFT;
F = CIF;
T = Tnominal + LGD;

% --------------- interesting stuff ends here ---------------

% -- assign return values, or make a pretty picture --
%
if nargout == 0

    % make a nice plot of the reassigned spectrum
    % using ratoimage, after cleaning up the data
    nhoriz = max( 500, size(X,2)*2 );
    nvert = max( 400, size(X,1)*2 );

    Fmax = .3*Fs; 
    Fmin = 0;
    Tmax = length(A) / Fs;
    Tmin = 0;
    ampthresh = -50;
    
    [X,F,T] = raframe( X, F, T, ampthresh, Fmin, Fmax, Tmin, Tmax );        
    ratoimage( X, F, T, nhoriz, nvert );

elseif nargout == 1,
    Bo = X;
elseif nargout == 2,
    Bo = X;
    Fo = F;
elseif nargout == 3,
    Bo = X;
    Fo = F;
    To = T;
end

% --------------- argument checking function ---------------

function [msg,x,nfft,Fs,win,nolap] = check_args(args)
%   Extract and check up to five arguments to findifspecgram
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
if (length(args) > 1) & ~isempty(args{2})
    nfft = args{2};
else
    nfft = min(length(x),256);
end
if (length(args) > 2) & ~isempty(args{3})
    Fs = args{3};
else
    Fs = 2;
end
if length(args) > 3 & ~isempty(args{4})
    win = args{4};
else
    if length(nfft) == 1
        win = hanning(nfft);
    else
        msg = 'findifspecgram: You must specify a window function.';
    end
end
if length(win) == 1
    win = hanning(win); 
end
if (length(args) > 4) & ~isempty(args{5})
    nolap = args{5};
else
    nolap = ceil(length(win)/2);
end
 
% do error checking
if (length(nfft)==1) & (nfft<length(win)),
    msg = 'findifspecgram: Window length cannot be greater than the FFT length.';
end
if (nolap >= length(win)),
    msg = 'findifspecgram: Overlap must be less than the window length.';
end
if (length(nfft)==1) & (nfft ~= abs(round(nfft)))
    msg = 'findifspecgram: FFT length must be a positive integer.';
end
if (nolap ~= abs(round(nolap))),
    msg = 'findifspecgram: Overlap must be a positive integer.';
end
if min(size(x))~=1,
    msg = 'findifspecgram: Input signal must be a vector (either row or column).';
end

% --------------- Sean Fulop's colormap ---------------

function J = myjet(m)
% This is adapted from jet.m by Sean Fulop.
% We use it to make nice-looking scatter plots
% of reassigned spectra.
if nargin < 1
   m = size(get(gcf,'colormap'),1);
end
n = ceil(m/4);
u = [(1:1:n)/n ones(1,n-1) (n:-1:1)/n]';
g = ceil(n/2) - (mod(m,4)==1) + (1:length(u))';
r = g + n;
b = g - n;
g(g>m) = [];
r(r>m) = [];
b(b<1) = [];
J = zeros(m,3);
J(r,1) = u(1:length(r));
J(g,2) = u(1:length(g));
J(b,3) = u(end-length(b)+1:end);

J(1,:) = [1,1,1];

