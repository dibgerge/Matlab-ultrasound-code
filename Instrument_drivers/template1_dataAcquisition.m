%% Start the generator and oscilloscope
import Lamby.*
fncgen = startGenerator();
oscilloscope = startOscilloscope();

%% Set the actuation signal
% --- Generate the excitation signal shape --- %
fc = 40e3;
BW = 0.5;
T = 10/fc;
[t, y, f, Y] = gauspulse(fc, BW, T, 2^14, 0, 0); 

% --- paramaters --- %
burst_period =  0.5; % seconds 
timebase = T/8;
time_off = 14*timebase;

% --- Setup the function generator --- %
set(fncgen, 'Frequency', 1/T);
set(fncgen, 'Amplitude', 10);
invoke(fncgen.Arbitrarywaveform, 'create', 'volatile', y);
invoke(fncgen.Arbitrarywaveform, 'setwaveform', 'volatile');
set(fncgen.burstmod, 'Cycles', 1)
set(fncgen.burstmod, 'InternalRate', 1/burst_period)
set(fncgen.burstmod, 'Enabled', 'on')
set(fncgen, 'Waveform', 'user');

%% Setup the oscilloscope 
set(oscilloscope.Acquire, 'AType', 'average')
% set(oscilloscope.Acquire, 'AType', 'normal')
set(oscilloscope.Acquire, 'Average', 16)
invoke(oscilloscope, 'Run');
set(oscilloscope.TriggerEdge, 'Source', 'ext');
set(oscilloscope.TriggerEdge, 'Sweep', 'auto');
set(oscilloscope.Waveform, 'Format', 'ASCII');
set(oscilloscope.Waveform, 'Mode', 'Raw');
set(oscilloscope.Waveform, 'Points', 20480);
set(oscilloscope.Timebase, 'Format', 'YT');
set(oscilloscope.Timebase, 'Offset',time_off)
set(oscilloscope.Timebase, 'Scale', timebase)

% --- Setup Channels--- %
set(oscilloscope.Channel1, 'BWLimit', 'on');
set(oscilloscope.Channel1, 'Probe', '1X');
set(oscilloscope.Channel1, 'Scale', 5);
set(oscilloscope.Channel2, 'BWLimit', 'on');
set(oscilloscope.Channel2, 'Probe', '10X');
set(oscilloscope.Channel2, 'Scale', 0.05);
% set(oscilloscope.Channel3, 'BWLimit', 'on');
% set(oscilloscope.Channel3, 'Probe', '1X');
% set(oscilloscope.Channel4, 'BWLimit', 'on');
% set(oscilloscope.Channel4, 'Probe', '10X');

%% --- Obtain data from oscilloscope --- %
nsamples = 10240;
invoke(oscilloscope, 'Stop');
t = zeros(nsamples,2);
meas = zeros(nsamples, 2);

for i=1:2
    set(oscilloscope.Waveform, 'Source', ['channel' num2str(i)]);
    dt = get(oscilloscope.Waveform, 'XIncrement');
    t(:,i) = 0:dt:(nsamples-1)*dt;
    dat = get(oscilloscope.Waveform, 'Data');
    meas(:,i) = sscanf(dat(12:end), '%f,');
end

plot(t, meas)

%% Disconnect
disconnect(oscilloscope);
disconnect(fncgen)

