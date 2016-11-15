function [idx] = findRidges( F, fs, nfft )
%   findRidges - return the indices corresponding to
%                the frequency ridges in a reassigned
%                spectral surface. Use this function
%                to find the spectral ridges corresponding
%                to concentrations of energy (like 
%                sinusoids) reassigned spectral data. 
% ---------------------------------------------------------
% [idx] = findRidges( F, Fs, NFFT )
%
% F =    the reassigned frequency matrix 
% Fs =   the sample rate
% NFFT = the size of the FFT used to compute the reassigned
%        spectrum
% 
% idx =  the indices of data in the reassigned matrices
%        corresponding to ridges oriented along the 
%        frequency dimension, that is, the indices 
%        corresponding to points of locally minimal 
%        frequency reassignment
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


[rows,cols] = size(F);
Fnom = ((0:rows-1)'*fs/nfft)*ones(1,cols);
fra = F - Fnom;
absfra = abs(fra);
ismin = absfra <= circshift( absfra, 1 ) & absfra <= circshift( absfra, -1 );
isswitch = circshift( fra, 1 ) > 0 & circshift( fra, -1 ) < 0;
idx = find( (ismin.*isswitch) == 1 );
