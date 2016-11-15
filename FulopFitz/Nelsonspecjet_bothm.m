
% Author: Sean Fulop

function [STFTpos, CIFpos, LGDpos, LGDderiv, f, t] = Nelsonspecjet_bothm(signal, Fs, window, overlap, fftn, low, high, clip)
  if nargin ~= 8
    error('usage: ([STFT [, CIF [, LGD [, f [, t]]]]] = Nelsonspecjet(signal [, Fs [, window [, overlap [, fftn [, low [, high [, clip]]]]]]]) ')
  end

  % Arguments: signal is the signal, with sampling rate Fs;
  % window is the analysis frame length, in sample points;
  % overlap is the number of points the frames overlap, equaling window minus frame
  % advance;
  % fftn is the discrete Fourier transform analysis length, in samples;
  % low and high are the frequency limits of the display, in Hz;
  % clip is a negative number which specifies the number of dB down from
  % the maximum that are to be shown---everything quieter is clipped

delay = 1;

len = length(signal);

freqdelay = 1;
% for now just use a Hanning window
window = hanning(window);


  
% compute window offsets
win_size = length(window);
step = win_size - overlap;

% fftn = max(512, win_size);
% sets fft frame length to at least 512; automatically zero filled;
% purpose of this is to achieve better frequency quantization so the Nelson
% spectrum has less sparse representation for each data window

% build matrices of delayed and undelayed windowed data slices
% each column is a data window
% Follow Paul Kienzle's Octave program here
offset = [ 0 : step : len-win_size-1 ];
WW = window*ones(1,length(offset));

%S = zeros (fftn, length(offset));
%Sdel = zeros (fftn, length(offset));
idx = ((1:win_size)'*ones(1,length(offset))) + (ones(win_size,1)*offset);
S = signal(idx+1).*WW;
Sdel = signal(idx).*WW;

% compute short-time fourier transform surfaces, having length(offset)
% columns
% each column of these is an fftn-length fft of a signal window
STFTdel = fft(Sdel,fftn);
STFT = fft(S,fftn);
STFTfreqdel = [ STFT(fftn,:); STFT(1:fftn-1,:) ];
%STFTfreqdel = [zeros(freqdelay, length(offset)); STFT(1:fftn-freqdelay, 1:end)];

% Double delayed surface used for cross-spectral estimates of mixed partial
% derivatives
%STFTfrtimedel = [zeros(freqdelay, length(offset)); STFTdel(1:fftn-freqdelay, 1:end)];  
STFTfrtimedel = [ STFTdel(fftn,:); STFTdel(1:fftn-1,:) ];

% extract the positive frequency components
if rem(fftn,2)==1
  ret_n = (fftn-1)/2;
else
  ret_n = fftn/2;
end

%STFTpos = STFT(1:ret_n, :);
%STFTdelpos = STFTdel(1:ret_n, :);
%STFTfreqdelpos = STFTfreqdel(1:ret_n, :);
%STFTfrtimedelpos = STFTfrtimedel(1:ret_n, :);

if high > Fs*(ret_n-1)/fftn
    high = Fs*(ret_n-1)/fftn;
end

lowindex = round(low/Fs*fftn);
if lowindex == 0
    lowindex = 1;
end
highindex = round(high/Fs*fftn);

% Truncate STFT matrices to the frequency range to be plotted

STFTpos = STFT(lowindex:highindex, :);
STFTdelpos = STFTdel(lowindex:highindex, :);
STFTfreqdelpos = STFTfreqdel(lowindex:highindex, :);
STFTfrtimedelpos = STFTfrtimedel(lowindex:highindex, :);

% compute Nelson's intermediate cross-spectral surfaces

C = STFTpos .* conj(STFTdelpos);

L = STFTpos .* conj(STFTfreqdelpos);

MixCIF = STFTpos .* conj(STFTdelpos) .* conj(STFTfreqdelpos .* conj(STFTfrtimedelpos));

MixLGD = STFTpos .* conj(STFTfreqdelpos) .* conj(STFTdelpos .* conj(STFTfrtimedelpos));

% compute channelized instaneous frequency, a matrix of remapped frequency vectors,
% one vector for each time step for which ffts were computed
argC = mod(angle(C),2*pi);

argMixCIF = mod(angle(MixCIF),2*pi);

argMixLGD = mod(angle(MixLGD),2*pi);

CIFpos = ((Fs/delay).* argC)./(2.*pi);

% compute local group delay, a matrix of time adjustments,
% one (row) vector for each discrete frequency step determined by fft length
%LGDpos = -((win_size/Fs) .* arg(L))./(2.*pi);
argL = mod(angle(L), -2*pi);

LGDpos = -((fftn/Fs) .* argL)./(2.*pi);

% compute mixed partial time derivative of CIF by Nelson method
CIFderiv = ((Fs/delay) .* argMixCIF .* (fftn/win_size) .* argMixCIF)./(2.*pi);

% compute mixed partial frequency derivative of LGD by Nelson method
LGDderiv = ((fftn/win_size) .* argMixLGD .* (Fs/delay) .* argMixLGD)./(2*pi);


%f = [0:ret_n-1]*Fs/fftn;
t = (offset + win_size/2) ./ Fs;
%t = offset ./ Fs;

for a = 1:highindex-lowindex+1
%    tremap(a,:) = LGDpos(a,:) + t + (fftn-win_size)./(Fs);
%    tremap(a,:) = LGDpos(a,:) + t + (fftn - (1/2)*win_size)./Fs;
   tremap(a,:) = LGDpos(a,:) + t - (win_size/2 - 1)./Fs;
end
 
STFTmag = abs(STFTpos); % magnitude of STFT
STFTmag = STFTmag ./ max(max(STFTmag)); % normalize so max magnitude will be 0 db
STFTmag = 20*log10(STFTmag); % convert to log magnitude, decibel scale

%clip = -30; % set dynamic range dB

plot_these = find(STFTmag >= clip & low <= CIFpos & CIFpos <= high & t(1) <= tremap & tremap <= t(end) & (abs(CIFderiv) < 0.25 | abs(LGDderiv - 1) < 0.25));
% Thanks to Kelly Fitz for this trick

            
STFTplot = STFTmag(plot_these);
CIFplot = CIFpos(plot_these);
tremap = tremap(plot_these);

f = Fs*[lowindex-1:highindex-1] ./ fftn;


%for x = 1:length(offset)
plotpoints = [1000*tremap(:), CIFplot(:), STFTplot(:)];
%plotall = [plotall,plotpoints];
%end

% create figure, specifying WindowButtonMotionFcn property to a function callback
% write function which displays frequency coord of figure property CurrentPoint
% in a separate text window
global measure2;

figure('Position',[400,300,900,400],'WindowButtonMotionFcn', @showcoords2);
%figure('Position',[400,300,900,400]);
%for x = 1:length(offset)
%scatter3(1000*tremap(:, x), CIFplot(:, x), STFTplot(:, x), 1, STFTplot(:, x), 'filled'), view(0, 90);
scatter3(plotpoints(:, 1), plotpoints(:, 2), plotpoints(:, 3), 2, plotpoints(:, 3), 'filled'), view(0, 90);

%axis xy;
%colormap(flipud(gray(256)));
colormap(myjet);
%colormap(flipud(autumn));

global specgm;
specgm = gca;
hold on;
%end

xlabel('Time (ms)','FontName','times','FontSize',14)
ylabel('Frequency (Hz)','FontName','times','FontSize',14)
set(gca,'FontName','FixedWidth','FontSize',12)
%Use following if manual time axis setting is desired
%set(gca,'XLimMode','manual')
%set(gca,'XLim',[0 (t(length(offset))+t(1))*1000])
%set(gca,'XLim',[35 45])
hold off;
global meas;
meas = figure('Position',[100,0,150,50],'MenuBar','none');
set(gca,'FontName','FixedWidth','FontSize',12,'Visible','off')
%set(gcf,'CurrentAxes',outer);
measure2 = text(0, .95, 'teststring');

inst = input('Do you want a regular spectrogram (y/n)?','s');
if strcmp(inst,'n') | strcmp(inst,'no')
    return
else

for slice = 1:length(offset)
    %for x = lowindex:highindex
    STFTmag(:,slice) = max(STFTmag(:,slice),clip); %clip spgm plot
        %end
end

tforspgm = (offset + win_size/2) ./ Fs;

figure('Position',[10,10,300,250]);
pcolor(1000*tforspgm, f, STFTmag );
axis xy;
colormap(flipud(gray));
shading interp;
xlabel('Time (ms)','FontName','times','FontSize',12)
ylabel('Frequency (Hz)','FontName','times','FontSize',12)
set(gca,'FontName','FixedWidth')

end

%global measure2;

function showcoords2(obj, eventdata)

%global outer;
%global inner;


pointerloc = get(specgm, 'CurrentPoint');
freq = num2str(pointerloc(1,2));
location = strcat('frequency: ', freq);
delete(measure2);
set(0,'CurrentFigure', meas);
measure2 = text(0, .95, location);
end

end