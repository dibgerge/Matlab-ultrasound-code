function [deriv] = hainsworthtime(varargin)
%   hainsworthtime -
%               compute the partial derivative of the short
%               time amplitude spectrum with respect to time,
%               as proposed by Stephen Hainsworth and Malcolm
%               Macleod. This is said to be equivalent, after
%               scaling, to time reassignment. See "Time-frequency
%               reassignment: a review and analysis", technical
%               report, Cambridge University Engineering Dept, 
%               2003.
%
%				This is pretty much untested.
%---------------------------------------------------------------------
%   deriv = hainsworthtime(x,NFFT,Fs,window,olap)
%
%   deriv =  Hainsworth and Macleod's time derivative of 
%            amplitude spectrum
%   x =      signal
%   NFFT =   number of points in spectrum computation
%   Fs =     sample rate
%   window = window function for spectrum computation
%   olap =   overlap of analysis windows (number of samples,
%            default half the window length)
%
%   See also hainsworthfreq.

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
%   Last updated 12 March 2007 by Fitz.
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

% --------------- interesting stuff begins here ---------------

% -- compute spectra --
X = specgram(A,NFFT,Fs,W,olap);
magsqrd = abs(X).^2;
nonzero = find(magsqrd>0);
Nw = length(W);

% -- construct reassignment windows --
if ( mod(Nw,2) )
    Mw = (Nw-1)/2;
    framp = [(0:Mw),(-Mw:-1)]';
    tramp = (-Mw:Mw)';
else
    Mw = Nw/2;
    framp = [(0:Mw-1),(-Mw:-1)]' + 0.5;
    tramp = (-Mw:Mw-1)' + 0.5;
end

% -- scale the ramps to the desired units --
tramp = tramp/Fs;       % ramp in seconds
framp = framp * Fs/Nw;  % ramp in Hz

Wdt = -imag(ifft(framp.*fft(W)));

% -- compute auxiliary spectrum --
Xdt = specgram(A,NFFT,Fs,Wdt,olap);

% -- compute the derivative of amplitude spectrum --
[rows,cols] = size(X);
deriv = zeros(rows,cols);
deriv(nonzero) = real(Xdt(nonzero).*conj(X(nonzero)))./magsqrd(nonzero); 

              
% --------------- argument checking function ---------------

function [msg,x,nfft,Fs,win,nolap] = check_args(args)
%   Extract and check up to five arguments to hainsworthtime
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
        msg = 'hainsworthtime: You must specify a window function.';
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
    msg = 'hainsworthtime: Window length cannot be greater than the FFT length.';
end
if (nolap >= length(win)),
    msg = 'hainsworthtime: Overlap must be less than the window length.';
end
if (length(nfft)==1) & (nfft ~= abs(round(nfft)))
    msg = 'hainsworthtime: FFT length must be a positive integer.';
end
if (nolap ~= abs(round(nolap))),
    msg = 'hainsworthtime: Overlap must be a positive integer.';
end
if min(size(x))~=1,
    msg = 'hainsworthtime: Input signal must be a vector (either row or column).';
end

