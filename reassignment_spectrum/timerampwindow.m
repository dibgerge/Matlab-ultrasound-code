function [Wr] = timerampwindow( W, fs, nramps )
% timerampwindow - compute a window function scaled by a time ramp 
% such as would be used for reassigned spectral analysis. The window
% function W, supplied, is scaled by a linear time ramp symmetric
% about its midpoint. This is equivalent to a window whose Fourier
% transform is the derivative of the Fourier transform of W.
%---------------------------------------------------------------------
% [Wr] = timerampwindow( W, fs, nramps )
%
% W =       an analysis window function, assumed to
%           be symmetrical about n=0
% fs =      the sampling rate (in Hz), used to convert the units
%           of time-frequency reassignment calculations made
%           with the frequency derivative window from samples
%           to seconds. The ramped window is scaled by 
%           1/fs^nramps. (default fs = 1)
% nramps =  the number of times to apply a time ramp,
%           (default nramps = 1)
% Wr =      the window function that is the result of
%           W multiplied by a time ramp nramps times
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
%   Last updated 23 Apr 2007 by Fitz.
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
    error( 'timerampwindow: the sample rate must be a positive scalar' )
end

if nargin < 3
    nramps = 1;
elseif nramps <= 0 || ~ isscalar(nramps)
    error( 'timerampwindow: the number of derivatives must be a positive scalar' )
end

% -- construct time ramp --
Nw = length(W);
if ( mod(Nw,2) )
    Mw = (Nw-1)/2;
    tramp = (-Mw:Mw)';
else
    Mw = Nw/2;
    tramp = (-Mw:Mw-1)' + 0.5;
end

% -- apply time ramp --
for i = 1:nramps
    W = tramp.*W;
end

% -- scale to correct units --
Wr = W / (fs^nramps);
