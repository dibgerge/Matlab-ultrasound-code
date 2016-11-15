
% computes the Nelson power spectrum and plots it together with the fft-computed
% power spectrum (log magnitude)
% precede use by
% [signal, Fs, bits] = wavread('yourwavfile.wav');
% this will read a wav file into a signal matrix, which Nelsonpower expects.


% The delay parameter is the number of time indices the signal should be delayed
% to compute the cross-spectrum with itself; it is currently fixed to 1.

% arguments to function are name of the signal matrix, its sample rate, the
% low and high limits to the frequency range for plotting, in Hz.

function [PS, CIF, f] = Nelsonpower(signal, Fs, low, high)

delay = 1;

len = length(signal);

fftn = 2^nextpow2(4.*(len+delay));
% sets length of fft at 4 times signal length; automatically zero-filled;
% purpose of this is to achieve better frequency quantization so the Nelson
% spectrum shows more points

W1 = hanning(len+delay).*[zeros(delay, 1); signal];
W2 = hanning(len+delay).*[signal; zeros(delay, 1)];
% produces two length(signal)+delay windowed signals by applying a length-point 
% hanning window to the signal, once delayed, once undelayed


C = fft(W2,fftn) .* conj(fft(W1,fftn));

% computes Nelson's intermediate cross-spectral surface,
% which he has called the "cross-spectrum",
% in this case just a complex vector because the time point is fixed.
% Nelson's epsilon is the delay parameter.
% the vector is a discrete function of (angular) frequency

CIF = ((Fs/delay).* arg(C))/(2.*pi);

% computes an estimate of the channelized instantaneous frequency at a fixed time point.
% the vector is a discrete function of (angular) frequency, and its values
% are frequencies in Hz.
% It estimates the instantaneous frequency of the signal at the fixed time point
% (really the end of the signal, we are not doing t-f analysis yet) for each
% frequency "channel".


  % extract the positive frequency components
  if rem(fftn,2)==1
    ret_n = (fftn-1)/2;
  else
    ret_n = fftn/2;
end

if high > Fs*(ret_n-1)/fftn
    high = Fs*(ret_n-1)/fftn;
end

PS = 10 .* log10(fft(W2,fftn) .* conj(fft(W2,fftn)));
%f = Fs*[0:ret_n-1]/fftn;
lowindex = round(low/Fs*fftn);
if lowindex == 0
    lowindex = 1;
end
highindex = round(high/Fs*fftn);
f = Fs*[lowindex:highindex]/fftn;
% computes the decibel-unit magnitude spectrum of the signal

for x = 1:ret_n
    if CIF(x) < low | CIF(x) > high
        CIFplot(x) = low;
        PSplot(x) = 0;
    else CIFplot(x) = CIF(x);
        PSplot(x) = PS(x);
    end
end

% create figure, specifying WindowButtonMotionFcn property to a function callback
% write function which displays x,y coords of figure property CurrentPoint
% in a separate text window
global measure;

global spectrum;
spectrum = figure('WindowButtonMotionFcn', @showcoords);
global outer;
global inner;
outer = axes('Position',[0 0 1 1],'Visible','off');
inner = axes('Position',[.1 .1 .7 .8]);
plot(CIFplot(1:ret_n), PSplot(1:ret_n), 'd', 'MarkerSize', 2, 'MarkerEdgeColor', 'r')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
%plot(CIFplot(1:ret_n), PSplot(1:ret_n), '-')
set(spectrum,'CurrentAxes',outer);
measure = text(0, .95, 'teststring');
%coords = figure


inst = input('Overlay the regular fft spectrum (y/n)?','s');
if strcmp(inst,'n') | strcmp(inst,'no')
    return
else
%    hold on
    set(spectrum,'CurrentAxes',inner);
    set(inner,'NextPlot','add');
    plot(f, PS(lowindex+1:highindex+1), 'g')
%    hold off
% plots the cross-spectrally computed Nelson spectrum on the same graph as the 
% standard log magnitude power spectrum
end

end

function showcoords(obj, eventdata)
global measure;
global outer;
global inner;
global spectrum;
pointerloc = get(inner, 'CurrentPoint');
freq = num2str(pointerloc(1,1));
amp = num2str(pointerloc(1,2));
location = strcat('frequency: ', freq, '  intensity: ', amp);
delete(measure);
set(spectrum,'CurrentAxes',outer);
measure = text(0, .95, location);
end
