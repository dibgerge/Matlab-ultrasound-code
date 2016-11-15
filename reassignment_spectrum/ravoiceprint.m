function [] = ravoiceprint( sig, fs, win, ampthresh, nhoriz, nvert )
% ravoiceprint - make and plot a reassigned voiceprint image
% --------------------------------------------------------------------
% ravoiceprint( sig, fs, win, threshdb, nhoriz, nvert )
%
% sig       = the signal
% fs        = the sample rate
% win       = the analysis window, or window length (in samples) to 
%             use for reassigned spectral analysis 
% threshdb  = the magnitude threshold (in dB), data points below
%             this threshold will not be plotted (default -50dB)
% nhoriz    = the number of pixels in the horizontal (time) dimension
%             (default 800)
% nvert     = the number of pixels in the vertical (frequency)
%             dimension (default 600)
%
% Returns nothing, just makes a pretty picture, suitable for framing.
%
% Uses: racomponents, raframe, ratoimage

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

% check arguments
if nargin < 3
    win = 0.005 * fs;
end
win = windowparam( win );

if nargin < 4
    ampthresh = -50;
end
if nargin < 5
    nhoriz = 800;
end
if nargin < 6
    nvert = 600;
end

% parameters
Fmax = 5000; 
Fmin = 0;
Tmax = length(sig) / fs;
Tmin = 0;
nfft = max( 256, 2 ^ ( ceil(log2( max( 2*length(win), nvert ) ) ) ) );
hop = max( 1, floor( length(sig) / nhoriz ) );
olap = max(1, length(win)-hop);

[X,F,T] = racomponents( sig, nfft, fs, win, olap, 0.2, 0.1 );
[X,F,T] = raframe( X, F, T, ampthresh, Fmin, Fmax, Tmin, Tmax );
ratoimage( X, F, T, nhoriz, nvert, ampthresh );


% --------------- window making function ---------------

function [W] = windowparam( arg )
% Return a short-time spectral analysis window from the specified
% argument. If the argument is a vector, return it as a
% column. If it is a scalar, return a Kaiser window of
% that length.
% -------------------------------------------------------

if length(arg) == 1
    W = kaiser(arg, 9); 
else
    W = arg;
end

W = W(:);

    
