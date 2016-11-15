function [XX,FF,TT] = raframe( X, F, T, xthreshdb, fmin, fmax, tmin, tmax )
% raframe - crop reassigned time-frequency data by applying 
%           an amplitude threshold, and time and frequency
%           boundaries. The data is NOT in matrix form when
%           it is returned.
% -----------------------------------------------------------------------
% [XX,FF,TT] = raframe( X, F, T, Xmindb, Fmin, Fmax, Tmin, Tmax )
%
%   X =      STFT matrix, columns are discrete Fourier 
%            transforms of short-time slices of the signal.
%   F =      frequency matrix containing reassigned 
%            frequencies (in Hz) for each point in X
%   T =      time matrix containing reassigned times
%            (in seconds) for each point in X
%   Xmindb = magnitude threshold (in dB), applied to X 
%   Fmin =   lower bound on frequency
%   Fmax =   upper bound on frequency
%   Tmin =   lower bound on time
%   Tmax =   upper bound on time
%
%   XX =     cropped STFT values
%   FF =     cropped freqencies corresponding to points in XX
%   TT =     cropped times corresponding to points in XX
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


Xmax = max(abs(X(:)));
M = 20*log10(abs(X)/Xmax);

inrange = find( F < fmax & F > fmin & T < tmax & T > tmin );
loud = find( M > xthreshdb );
select = intersect( inrange, loud );

XX = X(select);
FF = F(select);
TT = T(select);
