function [f, Y] = findfft(y, Fs, NFFT, option)
%FINDFFT computes the fourier transform and gives the frequency content
%Inputs:
%  y  Input signal, could be a vector or a matrix. fft is computed for each
%     column if it is a matrix.
%  Fs Sampling rate for y. If y is a matrix, the sampling rate for
%       all the columns should be the same.
%  NFFT number of points for the Fourier transform. If this is empty or not
%       provided, the length of y (or column length if y is a matrix)
%       determines the number of FFT points.
%  option: 'single', 'shift', or empty, corresponding to outputting the
%          sigle sided fft, or the shifted fft (true axes), or do nothing

if nargin > 2 && ~isempty(NFFT)
    if ~isnumeric(NFFT) || ~isscalar(NFFT) || NFFT <= 0 || mod(NFFT, 1) ~= 0
        error('NFFT should be a positive integer.')
    end
else
    NFFT = size(y,1);
end

df = Fs/NFFT;
f = (0:(NFFT-1))'*df;
f(f >= Fs/2) = f(f >= Fs/2) - Fs;
Y = fft(y, NFFT);
%f = fftshift(f);

if nargin == 4
    if strcmpi(option, 'shift')
        f = fftshift(f);
        Y = fftshift(Y);
    elseif strcmpi(option, 'single')
        f = f(1:ceil(end/2));
        if isvector(Y)
            Y = 2*Y(1:ceil(end/2))/length(Y);
        else
            Y = 2*Y(1:ceil(end/2),:)/size(Y,1);
        end
    end
end

end

