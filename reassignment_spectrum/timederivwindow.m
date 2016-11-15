function [Wd] = timederivwindow( W, fs, nderivs )
% timederivwindow - compute a window function that is the derivative
% of the supplied window, such as would be used for reassigned 
% spectral analysis. The derivative is computed by applying a frequency
% ramp to the Fourier transform of the supplied window function W.
%---------------------------------------------------------------------
% [Wd] = timederivwindow( W, fs, nderivs )
%
% W =       an analysis window function, assumed to
%           be symmetrical about n=0
% fs =      the sampling rate (in Hz), used to convert the units
%           of time-frequency reassignment calculations made
%           with the time derivative window from 1/samples
%           to 1/seconds (Hz). The derivative window is scaled
%           by fs^nderivs. (default fs = 1)
% nderivs = the number of time derivatives to compute,
%           (default 1)
% Wd =      the window function that is the result of
%           nderivs time derivatives of W
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
%   Last updated 23 April 2007 by Fitz.
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
if nargin < 2
    fs = 1;
elseif fs <= 0 || ~ isscalar(fs)
    error( 'timederivwindow: the sample rate must be a positive scalar' )
end

if nargin < 3
    nderivs = 1;
elseif nderivs <= 0 || ~ isscalar(nderivs)
    error( 'timederivwindow: the number of derivatives must be a positive scalar' )
end

% -- construct frequency ramp --
Nw = length(W);
if ( mod(Nw,2) )
    Mw = (Nw-1)/2;
    framp = [(0:Mw),(-Mw:-1)]';
else
    Mw = Nw/2;
    framp = [(0:Mw-1),(-Mw:-1)]' + 0.5;
end

framp = framp / Nw;

% -- apply frequency ramp --
for i = 1:nderivs
    W = -imag(ifft(framp.*fft(W)));
end

% -- scale to correct units --
Wd = W * (fs^nderivs);
