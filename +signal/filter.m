function yout = filter(y, Fs, cutoff, type)

if nargin < 4 
    type = 'lowpass';
end

[f, Y] = signal.findfft(y, Fs);
Y = signal.filterFrequencies(f, Y, cutoff, type);
yout = ifft(Y);

end

