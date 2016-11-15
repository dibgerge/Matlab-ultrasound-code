function [X,F,T] = crossspectrum(varargin)
%   crossspectrum - compute a reassigned spectrum
%   of a signal using Nelson's cross-spectral 
%   surfaces to approximate the partial derivatives.
% --------------------------------------------------------------
%   [X,F,T] = crossspectrum(A,NFFT,Fs,window)
%
%   Returns the complex reassigned spectrum of the signal
%   A in X, with a column of reassigned frequencies in F
%   and a column of reassigned times (offset from the 
%   center of the analysis window) in T. X,F, and T have
%   NFFT/2 rows for NFFT even, or (NFFT-1)/2 rows for NFFT 
%   odd (i.e. rows corresponding to frequencies below Fs/2).
%
%   A =      input signal (samples)
%   NFFT =   FFT size (default min of 256 and length of x)
%   Fs =     sample rate (default 2)
%   window = analysis window (samples, or hannning window 
%            length, default nfft)
%
%   X =      discrete Fourier transform of A.
%   F =      column of reassigned frequencies (in Hz)
%            for each point in X
%   T =      column of reassigned times (in seconds)
%            for each point in X
%
%   See also crossspecgram, raspecgram, raspectrum.

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
%   Last updated 19 Feb 2007 by Fitz.
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
error( nargchk(1,4,nargin) )
[msg,A,NFFT,Fs,W] = check_args( varargin );
error( msg )

A = A(:); % make it a column
W = W(:); % make it a column

% -- make sure that A is not too long, and invoke crossspecgram --
Nx = min(length(A), length(W)+1);
A = A(1:Nx);
[X,F,T] = crossspecgram(A,NFFT,Fs,W);

% -- time is referenced to center of window --
T = T - (0.5 * (length(W)-1)/Fs);

% --------------- argument checking function ---------------

function [msg,x,nfft,Fs,win] = check_args(args)
%   Extract and check up to four arguments to crossspectrum
%   from an array.
%
%   args =  array of arguments
%   x =     input signal (samples)
%   nfft =  FFT size (default min of 256 and length of x)
%   Fs =    sample rate (default 2)
%   win =   analysis window (samples, or hannning window 
%           length, default nfft)
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
        msg = 'crossspectrum: You must specify a window function.';
    end
end
if length(win) == 1
    win = hanning(win); 
end
 
% do error checking
if (length(nfft)==1) & (nfft<length(win)),
    msg = 'crossspectrum: Window length cannot be greater than the FFT length.';
end
if (length(nfft)==1) & (nfft ~= abs(round(nfft)))
    msg = 'crossspectrum: FFT length must be a positive integer.';
end
if min(size(x))~=1,
    msg = 'crossspectrum: Input signal must be a vector (either row or column).';
end
