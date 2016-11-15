function [X,F,T] = racomponents(varargin)
%   racomponents - return raspecgram pruned using mixed 
%                  partial derivative of phase to find
%                  sines and impulses, within specified
%                  tolerances. Computes the reassigned
%                  spectrogram and prunes the data. Since
%                  only sinusoids and impulses are retained,
%                  the data is NOT in matrix form when it 
%                  is returned.
%---------------------------------------------------------------------
%   [X,F,T] = racomponents(A,NFFT,Fs,window,olap,sinetol,imptol,thresh)
%
%   A =      input signal (samples)
%   NFFT =   FFT size (default min of 256 and length of A)
%   Fs =     sample rate (default 2)
%   window = analysis window (samples, or hannning window 
%            length, default nfft)
%   olap =   overlap of analysis windows (number of samples,
%            default half the window length)
%   sinetol =tolerance for detecting sinusoids (deviation of 
%            mixed partial derivative of phase from -1)
%   imptol = tolerance for detecting impulses (deviation of 
%            mixed partial derivative of phase from 0)
%   thresh = amplitude threshold (in dB), components of 
%            magnitude less than thresh are ignored
%
%   X =      STFT data at selected time-frequency coordinates 
%            satisfing the component selection criteria
%            (the tolerances and the magnitude threshold).
%   F =      frequency matrix containing reassigned 
%            frequencies (same units as Fs) for each point in X
%   T =      time matrix containing reassigned times
%            (same units as 1/Fs) for each point in X

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
error( nargchk(1,8,nargin) )
[msg,A,NFFT,Fs,W,olap,sintol,imptol,ampthreshdb] = check_args( varargin );
error( msg )

% -- reassigned spectral analysis --
%[X,F,T] = raspecgram( A, NFFT, Fs, W, olap );
[mpd,X,F,T] = mixedpd( A, NFFT, Fs, W, olap );
Xmax = max(abs(X(:)));
M = 20*log10(abs(X)/Xmax);

% prune data
loud = find( M > ampthreshdb );
nelson = find(abs(mpd) < imptol | abs(1+mpd) < sintol);

show = intersect( nelson, loud );
 
X = X(show);
F = F(show);
T = T(show);


% --------------- argument checking function ---------------

function [msg,x,nfft,Fs,win,nolap,sintol,imptol,thresh] = check_args(args)
%   Extract and check up to five arguments to racomponents
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
if (length(args) > 1) && ~isempty(args{2})
    nfft = args{2};
else
    nfft = min(length(x),256);
end
if (length(args) > 2) && ~isempty(args{3})
    Fs = args{3};
else
    Fs = 2;
end
if length(args) > 3 && ~isempty(args{4})
    win = args{4};
else
    if length(nfft) == 1
        win = hanning(nfft);
    else
        msg = 'mixedpd: You must specify a window function.';
    end
end
if length(win) == 1
    win = hanning(win); 
end
if (length(args) > 4) && ~isempty(args{5})
    nolap = args{5};
else
    nolap = ceil(length(win)/2);
end
if (length(args) > 5) && ~isempty(args{6})
    sintol = args{6};
else
    sintol = 0.1;
end
if (length(args) > 6) && ~isempty(args{7})
    imptol = args{7};
else
    imptol = 0.1;
end
if (length(args) > 7) && ~isempty(args{8})
    thresh = args{8};
else
    thresh = -120;
end
 
% do error checking
if (length(nfft)==1) && (nfft<length(win)),
    msg = 'mixedpd: Window length cannot be greater than the FFT length.';
end
if (nolap >= length(win)),
    msg = 'mixedpd: Overlap must be less than the window length.';
end
if (length(nfft)==1) && (nfft ~= abs(round(nfft)))
    msg = 'mixedpd: FFT length must be a positive integer.';
end
if (nolap ~= abs(round(nolap))),
    msg = 'mixedpd: Overlap must be a positive integer.';
end
if min(size(x))~=1,
    msg = 'mixedpd: Input signal must be a vector (either row or column).';
end
