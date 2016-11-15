classdef utacquisition < handle
    %UTACQUISITION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        step_size = [];
        ref_channel = 1;
        data_channel = 2;
        s;
        t;
    end
    
    properties(Dependent=true, SetAccess=private, GetAccess=public)
        Ts;
        nsig;
    end
    
    methods
        % constructor
        function this = utacquistion()    
        end
        
        %__________________________________________________________________%
        % read files from lecroy oscilloscope
        function readfromlecroy(this, pth, fname, nmax)
            files = dir(pth);
            nm = {files.name};            

            signame = strfind(nm, ['C' num2str(this.data_channel) fname]);
            nsig = length([signame{:}]);

            if ~isempty(this.s)
                warning('The previous signal array will be overwritten.');
            end
            
            if nargin == 4
                if nsig > nmax
                    nsig = nmax;
                end
            end
            this.s = cell(nsig, 1);

            for i=1:nsig
                allname = sprintf(['C' num2str(this.data_channel) fname '%05d.trc'], i-1);
                wave = ReadLeCroyBinaryWaveform([pth allname]);
                this.s{i} = UTLib.utsignal(wave.x, wave.y);
                
                if ~isempty(this.ref_channel) && this.ref_channel ~= 0
                    allname = sprintf(['C' num2str(this.ref_channel) fname '%05d.trc'], i-1);
                    wave = ReadLeCroyBinaryWaveform([pth allname]);
                    this.s{i}.yref = wave.y;
                end
            end       

            if nsig == 0 
                error('No files found.');
            else
                this.resample(this.s{1}.Ts);
            end
        end

        
        %__________________________________________________________________% 
        function Y = fft2(this)
            y = zeros(length(this.t), length(this.s));
            for i=1:length(this.s)
                y(:,i) = this.s{i}.y;
            end
            
            Y = fft2(y);
%             Y = Y(1:ceil(end/2),1:ceil(end/2));
        end

        %__________________________________________________________________%        
        function resample(this, rate)            
            for i=1:length(this.s)
                this.s{i}.resample(rate);
            end
            this.t = this.s{1}.t;  
        end
        
        
        %__________________________________________________________________%
        function truncate(this, time)
            for i=1:length(this.s)
                this.s{i}.truncate(time);
            end
        end
                
        
        %__________________________________________________________________%
        % -- remove the mean of the signal
        function removemean(this)
            for i=1:length(this.s)
                this.s{i}.removemean();
            end
        end
        
        %__________________________________________________________________%   
        % -- filter the signal using ideal filter
        function filter(this, cutoff, type)
            for i=1:length(this.s)
                this.s{i}.filter(cutoff, type);
            end
        end
        
        %__________________________________________________________________%
        % -- find the spectrogram of the signals
        function [f, t, p] = spectrogram(this, window, overlap, nfft)
            p = cell(this.nsig, 1);
            for i=1:this.nsig
                [~, f, t, p{i}] = spectrogram(this.s{i}.y, window, overlap, nfft, 1/this.Ts);
            end
        end        
        
        
        %__________________________________________________________________%
        function zcross = tof(this, varargin)
            threshold = -12;
            udb = 20*log10(abs(this.s{1}.y/max(abs(this.s{1}.y))));
            ind = find(udb > threshold, 1);
            tstart = this.t(ind);
            clr = 'mgbkyc';
            for i=1:this.nsig
                x = diff(sign(this.s{i}.y));
                indx_up = find(x > 0);              
                indi = find(this.t(indx_up) > tstart, 2);
                zcross(i,:) = this.t(indx_up(indi));
                tstart = zcross(i,1);
%                 hold on, plot(this.t, this.s{i}.y, ['-' clr(1+mod(i-1, length(clr)))])
%                 plot(zcross(i), this.s{i}.y(indx_up(indi)), 'or');                
%                 xlim([15e-6 90e-6])
            end
        end
        
        %__________________________________________________________________%
        function times = dispersion(this, varargin)
            %% let the window be half the inverse of bandwidth
            [B, fmin, fmax] = SignalProc.findbandwidth(this.t, this.s{1}.y, -6);
            nsamples = floor((0.5/B)/this.Ts);
            nolap = floor(nsamples - 0.025*nsamples);
            
            if nsamples < 2^14
                nfft = 2^14;
            else
                nfft = nsamples;
            end            
            
            for i=1:length(this.s)
                [~, f, t, p] = spectrogram(this.s{i}.y, nsamples, nolap, nfft, 1/this.Ts);                           
                p = p/max(abs(p(:)));            
                ind1 = find(f > fmin/4, 1);
                ind2 = find(f > fmax*2, 1);
                fi = f(ind1:ind2);
                pi = p(ind1:ind2,:);
                if i==1
                    h = figure('Renderer', 'zbuffer', 'position', get(0,'Screensize'));
                    surf(t*1e6, fi*1e-3, 10*log10(pi)), view(2)
                    shading interp
                    ylim([fmin*1e-3/4 2*fmax*1e-3])
                    xlim([10 300])
                    caxis([-40 0])
                    grid off
                    xlabel('Time (\mu s)')
                    ylabel('Frequency (kHz)');
                    rect = getrect;
                    
                    indf1 = find(f > rect(2)*1e3, 1);
                    indf2 = find(f > 1e3*(rect(2)+rect(4)), 1);
                    indt1 = find(t > rect(1)*1e-6, 1);
                    indt2 = find(t > 1e-6*(rect(1)+rect(3)), 1);
                    winp = p(indf1:indf2,indt1:indt2);
                    [maxc, ~] = max(winp);
                    [~, Ir] = max(maxc);
                    close(h);
                    h = figure('Renderer', 'zbuffer', 'position', get(0,'Screensize'));
                    surf(10*log10(winp)), view(2), shading interp
                else             
                    winp = p(indf1:indf2,indt1:indt2);
                    [maxic, ~] = max(winp);
                    [~, Iir] = max(maxic);
                    nshift = Iir - Ir;
                    winp = p(indf1:indf2,(nshift+indt1):(nshift+indt2)); 

                    h = figure('Renderer', 'zbuffer', 'position', get(0,'Screensize'));
                    surf(10*log10(winp)), view(2), shading interp                
                end
            end
        end
        
        
        %__________________________________________________________________%
        function plot(this, num, varargin)
            if num == 0 
                num = 1:length(this.s);
            end
            
            for i=num
                this.s{i}.plot(varargin{:})
                hold all                
            end
        end
        
        %__________________________________________________________________%
        function Ts = get.Ts(this)
            Ts = mean(diff(this.t));
        end      
        
        function nsig = get.nsig(this)
            nsig = length(this.s);
        end              
        
    end
    
end

