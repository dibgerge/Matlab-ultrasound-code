function [Zo,Fo,To] = ratoimage( X, F, T, nhoriz, nvert, ampthresh )
% ratoimage - build an image (two-dimensional matrix)
%             depicting a reassigned spectrum by sampling
%             the reassigned data at regular intervals
%             in frequency and time. 
%---------------------------------------------------------------------
% [Z,Fk,Tn] = ratoimage( X, F, T, nhoriz, nvert [, ampthresh] )
%
% X 	= the reassigned spectral data
% F 	= the matrix of reassigned frequencies 
%         (same size as X)
% T 	= the matrix of reassigned times
%         (same size as X)
% nhoriz = the number of pixels in the horizontal (time) dimension
% nvert = the number of pixels in the vertical (frequency)
%
% ampthresh = an optional threshold (in dB) on the values of X, 
%             (absolute) values below this threshold are colored
%             like the background, if not specified no threshold
%             is applied
%
% Z 	= the sampled spectral data
% Fk 	= the sampled frequencies (same size as Z)
% Tn 	= the sampled times (same size as Z)
%
% If no output arguments are provided, the image is rendered
% as a log-magnitude plot (with a threshold 80 dB below the 
% maximum absolute value of X, if no threshold was specified).

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

if nhoriz < 1 | nvert < 1
    error 'horizontal and vertical image dimensions must be > 1'
end

T = T(:);
F = F(:);
X = X(:);

Tmax = max( T );
Tmin = min( T );
dt = ( Tmax - Tmin ) / ( nhoriz - 2 );
nmax = ceil( Tmax / dt );
nmin = floor( Tmin / dt );

Tn = Tmin + ( dt * (0:nmax - nmin) );

Fmax = max( F );
Fmin = min( F );
df = ( Fmax - Fmin ) / ( nvert - 2 );
kmax = ceil( Fmax / df );
kmin = floor( Fmin / df );

Fk = Fmin + ( df * (0:kmax - kmin) );

Z = zeros( nvert, nhoriz );

for i = 1:length(X)
	n = 1 - nmin + ( T(i) / dt );
	k = 1 - kmin + ( F(i) / df );
	alpha = n - floor(n);
	beta = k - floor(k);
	Z(floor(k),floor(n)) = Z(floor(k),floor(n)) + ((1-alpha) * (1-beta) * X(i));
	Z(ceil(k),floor(n)) = Z(ceil(k),floor(n)) + ((1-alpha) * (beta) * X(i));
	Z(floor(k),ceil(n)) = Z(floor(k),ceil(n)) + ((alpha) * (1-beta) * X(i));
	Z(ceil(k),ceil(n)) = Z(ceil(k),ceil(n)) + ((alpha) * (beta) * X(i));
end

% apply an amplitude threshold to the data
if nargin > 5
    zmin = 10^(0.05*(ampthresh));
    zbg = 10^(0.05*(ampthresh-10)); 
        % background level, 10dB below threshold
        
    % make values below the threshold the color of
    % the background
    Z( find( abs(Z) <= zmin ) ) = zbg;
end

% if no output arguments, display the image
if nargout == 0
    
    if nargin < 6
        % display 80 dB of data
        ampthresh = 20*log10(max(abs(Z(:)))) - 80;
    
        zmin = 10^(0.05*(ampthresh));
        zbg = 10^(0.05*(ampthresh-10)); 
            % background level, 10dB below threshold
        Z( find( abs(Z) <= zmin ) ) = zbg;
    end

    % display the image
    imagesc(Tn, Fk, 20*log10(abs(Z)));
    axis([min(Tn),max(Tn),min(Fk),max(Fk)]); 
    axis xy;
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    %colormap( 1 - gray )
        % grayscale looks better than the default jet
   
else
    Zo = Z;
    To = Tn;
    Fo = Fk;
end
