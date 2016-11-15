function [ deriv ] = dphasedf2(varargin)
%   dphasedf2 - compute the second partial derivative of 
%               spectral phase w.r.t. frequency
% ---------------------------------------------------------
%   deriv = dphasedf2(A,NFFT,Fs,window,olap)
%
%   deriv =  second derivative of phase w.r.t frequency
%   A =      signal
%   NFFT =   number of points in spectrum computation
%   Fs =     sample rate
%   window = window function for spectrum computation
%   olap =   overlap of analysis windows (number of samples,
%            default half the window length)
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

% -- check arguments --
% use the auxiliary function, below, 
error( nargchk(1,5,nargin) )
[msg,A,NFFT,Fs,W,olap] = check_args( varargin );
error( msg )

% -- compute spectrum --
Nx = min(length(A), length(W));
Nw = length(W);
[X,F,T] = specgram(A,NFFT,Fs,W,olap);
magsqrd = abs(X).^2; % used to compute time and freq corrections

% construct reassignment windows
Wt = timerampwindow(W,Fs);
Wtt = timerampwindow(Wt,Fs);

% -- compute auxiliary spectra --
Xt = specgram(A,NFFT,Fs,Wt,olap);
Xtt = specgram(A,NFFT,Fs,Wtt,olap);

% -- compute derivative --
magsqrd = max( magsqrd, eps );  % make sure that we don't divide by zero
    % this tiny adjustment can never be significant,
    % since the only way magsqrd could be 0 is if
    % X is zero, and in that case, conj(X), and hence
    % the numerators in the expressions below, will also
    % be zero, so this is just making sure that 0/0 gives 0.

term1 = imag( ( Xt .* conj( X ) ) .^ 2 ) ./ ( magsqrd .^ 2 );
term2 = imag( Xtt .* conj( X ) ) ./ magsqrd;
deriv = term1 - term2;


% --------------- argument checking function ---------------

function [msg,x,nfft,Fs,win,nolap] = check_args(args)
%   Extract and check up to five arguments to dphasedf2
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
        msg = 'dphasedf2: You must specify a window function.';
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
    msg = 'dphasedf2: Window length cannot be greater than the FFT length.';
end
if (nolap >= length(win)),
    msg = 'dphasedf2: Overlap must be less than the window length.';
end
if (length(nfft)==1) & (nfft ~= abs(round(nfft)))
    msg = 'dphasedf2: FFT length must be a positive integer.';
end
if (nolap ~= abs(round(nolap))),
    msg = 'dphasedf2: Overlap must be a positive integer.';
end
if min(size(x))~=1,
    msg = 'dphasedf2: Input signal must be a vector (either row or column).';
end

