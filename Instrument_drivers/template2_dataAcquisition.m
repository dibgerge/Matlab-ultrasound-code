%% Start the generator and oscilloscope
import Lamby.*
tic
fncgen = startGenerator();
oscilloscope = startOscilloscope();
date = '2014_07_08';
fall = [25e3 40e3 60e3 100e3];

fname = 'actuator_2B_sens_6B_7B_8B';
impact = '7';


nchannel = 4;
nsamples = 10240;

%% Setup the oscilloscope 
%set(oscilloscope.Acquire, 'AType', 'average')
set(oscilloscope.Acquire, 'AType', 'normal')
%set(oscilloscope.Acquire, 'Average', 2)
set(oscilloscope.TriggerEdge, 'Source', 'ext');
set(oscilloscope.TriggerEdge, 'Sweep', 'auto');
set(oscilloscope.Waveform, 'Format', 'ASCII');
set(oscilloscope.Waveform, 'Mode', 'Raw');
set(oscilloscope.Waveform, 'Points', 20480);
set(oscilloscope.Timebase, 'Format', 'YT');

% --- Setup Channels--- %
set(oscilloscope.Channel1, 'BWLimit', 'on');
set(oscilloscope.Channel1, 'Probe', '1X');
set(oscilloscope.Channel2, 'BWLimit', 'on');
set(oscilloscope.Channel2, 'Probe', '10X');
set(oscilloscope.Channel3, 'BWLimit', 'on');
set(oscilloscope.Channel3, 'Probe', '10X');
set(oscilloscope.Channel4, 'BWLimit', 'on');
set(oscilloscope.Channel4, 'Probe', '10X');

 burst_period =  0.2; % seconds 
 BW = 0.5;
 set(fncgen.burstmod, 'Cycles', 1)
 set(fncgen.burstmod, 'InternalRate', 1/burst_period)
 set(fncgen.burstmod, 'Enabled', 'on')
 set(fncgen, 'Waveform', 'user');
 set(fncgen, 'Amplitude', 10);
 
set(oscilloscope.Channel1, 'Scale', 5);
set(oscilloscope.Channel2, 'Scale', 0.15);
set(oscilloscope.Channel3, 'Scale', 0.15);
set(oscilloscope.Channel4, 'Scale', 0.15);

%% Begin data acquisition
for j=1:length(fall)
    fc = fall(j);
   
    T = 10/fc;
    [t, y, f, Y] = gauspulse(fc, BW, T, 2^14, 0, 0); 

    % --- paramaters --- %   
    timebase = T/8;
    time_off = 10*timebase;
    invoke(fncgen.Arbitrarywaveform, 'create', 'volatile', y);
    invoke(fncgen.Arbitrarywaveform, 'setwaveform', 'volatile');
    set(fncgen, 'Frequency', 1/T);
  
    % --- Obtain data from oscilloscope --- %   
    invoke(oscilloscope, 'Run');
    set(oscilloscope.Timebase, 'Offset',time_off)
    set(oscilloscope.Timebase, 'Scale', timebase)
    pause(1);
   
    invoke(oscilloscope, 'Stop');
    meas = zeros(nsamples, nchannel);
    
    for i=1:nchannel
        set(oscilloscope.Waveform, 'Source', ['channel' num2str(i)]);
        dt = get(oscilloscope.Waveform, 'XIncrement');
        t = 0:dt:(nsamples-1)*dt;       
        dat = get(oscilloscope.Waveform, 'Data');
        meas(:,i) = sscanf(dat(12:end), '%f,');
    end
    
   % figure;plot(t, meas); title([num2str(fc/1e3) 'khz'])
    
     save(['Data/' date '/' fname '_f' num2str(fc/1e3) 'khz_impact' impact], 't', 'meas');
end

%% Disconnect
disconnect(oscilloscope);
disconnect(fncgen)
toc
