classdef utsignal < UTLib.utmeasurement
    %UTSIGNAL provides operations needed for computing parameters from
    %ultrasonic signals
    % inherits from class utmeasurement
    
    properties(Access=public)
        yref = [];      % optional reference signal (could be the excitation pulse)
        y = [];         % the main UT time domin signal
        t = [];         % the time vector corresponding to the signal        
    end
    
    properties(Dependent=true, SetAccess=private, GetAccess=public)
        yenv;           % signal envelop
        ydb;            % signal in dB
        ny;             % normalized signal
        yenvdb;         % signal envelop in dB
        nyenv;          % normalized signal envelop 
        nydb;           % normalized signal in dB
        nyenvdb;        % normalized signal envelop in dB
        Ts;             % sampling rate
    end
    
    properties(Access=protected)
        options;
    end
    
    methods
        %__________________________________________________________________%
        function this = utsignal(t, y, varargin)            
            this@UTLib.utmeasurement(varargin{:});            
            this.y = y(:);            
            if ((length(t) ~= length(y)) && ~isscalar(t))
                error(['t and y should have the same length, or t ' ...
                'should be a scalar representing the sampling rate.']);            
            end
            
            if isscalar(t)
                this.t = (0:length(y)-1)*t;
            else
                this.t = t(:);
            end      
            
            % set the default threshold to be 6 dB from the mean of the signal
            thresh = mean(this.nyenvdb) + 13;
            % - method : 'amplitude', 'envelop', or 'correlation'          
            %            used for computing the time of flight
            this.options = struct('threshold', thresh, ...
                         'holdoff', 0, ...
                         'width', 2e-6, ...
                         'method', 'correlation');  
        end
        
        
        %__________________________________________________________________%
        function resample(this, ts)            
            t2 = this.t(1):ts:this.t(end);            
            this.y = interp1(this.t, this.y, t2);                            
                        
            if ~isempty(this.yref)
                this.yref = interp1(this.t, this.yref, t2);
            end
            this.t = t2;
        end
        
        
        %__________________________________________________________________%
        function [f, Y] = fft(this)
            [f, Y] = SignalProc.findfft(this.t, this.y, [], 'shift');
        end
        
        
        %__________________________________________________________________%
        % -- remove the initial bang from the signal
        function removeInitBang(this)
            thresh = mean(this.nyenvdb);
            maxloc = find(this.nyenvdb > 0, 1);
            endloc = find(this.nyenvdb(maxloc:end) < thresh, 1);
            this.y(1:(maxloc+endloc)) = this.y((maxloc+endloc));
        end
        
        
        %__________________________________________________________________%
        % -- filter the signal using ideal filter
        function filter(this, cutoff, type)
            [this.t, this.y] = SignalProc.filterSignal(this.t, this.y, cutoff, length(this.y), type);
        end
        
        
        %__________________________________________________________________%
        % -- remove the mean of the signal
        function removemean(this)
            this.y = this.y - mean(this.y);
        end
        
        
        %__________________________________________________________________%
        function truncate(this, time)
            ind = find(this.t > time);
            if ~isempty(ind)
                this.t = this.t(1:ind);
                this.y = this.y(1:ind);
            else
                warning('Given time is longer than signal.')
            end
        end
        
        
        %__________________________________________________________________%
        % -- get the time of flight the ultrasound signal ---
        function varargout = gettof(this,  varargin)
            % set given options
            this.setOptions(varargin);
            
            yseg = this.segment(varargin{:}); %#ok<*PROP>
            N = size(yseg, 2);
            tof = zeros(N-1,1);
            if strcmpi(this.options.method, 'amplitude')                 
                [~, I] = max(abs(yseg));                
                tof = diff(this.t(I));                    
            elseif strcmpi(this.options.method, 'envelop')
                [~, I] = max(abs(hilbert(yseg)));
                tof = diff(this.t(I));
            elseif strcmpi(this.options.method, 'correlation')                
                for i=2:N
                    packet1 = yseg(:,i-1)/max(abs(yseg(:,i-1)));
                    packet2 = yseg(:,i)/max(abs(yseg(:,i)));
                    c = xcorr(packet1, packet2, 'coeff');
                    [~,I] = max(c);
                    mx = length(packet1)-I;
                    tof(i-1) = mx*this.Ts;
                end
            end
            if nargout == 1
                varargout{1} = tof;
            elseif nargout == 2
                varargout{1} = tof;
                varargout{2} = yseg;
            end
        end
        
        
        %__________________________________________________________________%
        function [att, rmse] =  getattenuation(this, varargin)
            ynew = this.diffractioncorrect(varargin{:});            
            [mx, ind] = max(ynew);
            Dmx  = this.t(ind)*this.velocity;
            [f, gof] = fit(Dmx(:), mx(:), 'exp1');
            att = -8.686*f.b/100;
            rmse = gof.rmse;            
        end
        
        
        %__________________________________________________________________%
        function ynew = diffractioncorrect(this, varargin)
            [tof, yseg] = this.gettof(varargin{:});            
            Nseg = size(yseg, 2);
            % get the --approximate-- tof of the first wave packet (reference)
            [~,ind] = max(yseg(:,1));
            reftof = this.t(ind);
            
            % assume all values in 'tof' are close, small error. find the mean.
            meantof = mean(tof);            
            % using mean tof, find the number of reflections of the reference signal. 
            tstart = this.gettriggertime();
            N1 = round((reftof-tstart)/meantof);
            
            % calculate the mean velocity
            this.velocity = mean(2*this.thickness./tof);
            
            D = (N1:(N1+Nseg-1))*this.thickness*2; 
            ynew = zeros(size(yseg));
            for i=1:Nseg
                [this.frequency, Y] = SignalProc.findfft(this.t, yseg(:,i));                
                Dhat = this.diffraction(D(i));%exp(-sqrt(-1)*D(i)*this.k).*
                yrat = SignalProc.wdeconv(Y, Dhat, max(abs(Y))/10e4);
                ynew(:,i) = real(ifft(yrat));                
            end            
        end
        
        
        %__________________________________________________________________%
        % -- cut the UT signal into its component reflections ---
        function yseg = segment(this, varargin)         
            this.setOptions(varargin);

            pkloc = UTLib.peakfinder(this.nyenvdb, 6, this.options.threshold, 1, false);            
            % remove peaks before the holdoff time            
            pktime = this.t(pkloc);
            ind = find(pktime < (this.options.holdoff - this.t(1)));
            pkloc(ind) = [];    
            pktime(ind) = [];   

            % remove points within the same wave packet
            inprox = UTLib.proximityCheck(pktime, this.options.width/2);
            for i=1:length(inprox)
                % keep the maximum of all points in proximity
                [~, I] = max(this.nyenvdb(pkloc(inprox{i}(1):inprox{i}(end))));
                inprox{i}(I) = [];
            end            
            pkloc([inprox{:}]) = [];

            % remove end points if smaller than width 
            N = ceil(this.options.width/(2*this.Ts));            
            if pkloc(1) < N
                pkloc(1) = [];
            end
            if pkloc(end)+N > length(this.y)
                pkloc(end) = [];
            end            

            yseg = repmat(this.y, 1, length(pkloc));

            for i=1:length(pkloc)            
                if pkloc(i)-N > 0       % shouldnt happen, but keep for safety
                    yseg(1:pkloc(i)-N, i) = 0;%yseg(pkloc(i)-N);
                end
                if pkloc(i)+N < length(this.y)% shouldnt happen, but keep for safety
                    yseg(pkloc(i)+N+1:end, i) =0;% yseg(pkloc(i)+N+1);            
                end                                    
%                 wind = UTLib.shiftedHamming(length(this.y), pkloc(i));                   
%                 yseg(:,i) = yseg(:,i).*wind(:);
            end         
        end
        
        
        %__________________________________________________________________%        
        function tstart = gettriggertime(this)
            % if yref is given, get the exact triggering time
            % otherwise, assume t=0 is the trigger time.
            if ~isempty(this.yexc)
                nyref = abs(this.yexc)/max(abs(this.yexc));
                mref = 20*log10(mean(nyref));
                crossInd = find(20*log10(nyref) > mref+6, 1);   % index of  threshold crossing.
                tstart = this.texc(crossInd);
            else
                tstart = 0;
            end            
        end
        
        
        %__________________________________________________________________%
        function [pkind, pkmag] = getpeaks(this, varargin)            
            options = struct('sel', (max(this.yenv)-min(this.yenv))/4, ...
                'thresh', 0, ...
                'extrema', 1, ...
                'include_endpoints', 1);
            optionNames = fieldnames(options);
            nArgs = length(varargin);
            if round(nArgs/2)~=nArgs/2
                error('propertyName/propertyValue pairs are required.')
            end
            for pair = reshape(varargin,2,[])
                if any(strncmp(pair{1}, optionNames, length(pair{1})))
                    options.(pair{1}) = pair{2};
                else
                    error('%s is not a recognized parameter name.', pair{1})
                end
            end
            [pkind, pkmag] = UTLib.peakfinder(this.yenv, options.sel, options.thresh, ...
                options.extrema, options.include_endpoints);
        end
        
        
        %__________________________________________________________________%
        function [t, y] = deconvolve(this, varargin)
            options = struct('Type', 'Gaussian', ...
                        't', [], ...
                        'signal', [], ...
                        'Frequency', 5e6, ...
                        'Bandwidth', 0.28758);
            this.checkoptions(options, varargin);
            
            if strcmp(options.Type, 'Gaussian') && isempty(options.t) && isempty(options.signal)
                [tdeconv, ydeconv, ~, ~] = SignalProc.gauspulse(options.Frequency, ...
                    options.Bandwidth, 'T', this.t(end));
            end
            
            [t, yr, yne] = UTLib.syncSignals(this.t, this.y, tdeconv, ydeconv);
            [fr, Yr] = SignalProc.findfft(t, yr);            
            Yr = Yr/length(Yr);
            [fne, Yne] = SignalProc.findfft(t, yne);
            Yne = Yne/length(Yne);
            y = ifft(Yr.*Yne)*100000000*length(Yne);
        end
               
        
        %__________________________________________________________________%  
        function plot(this, varargin)            
            options = struct('xscale', 1e-6, ...
                'yscale', 1e-3, ...
                'yunit', 'Volts', ...
                'which', 1, ...     % two bit flag: 1) signal, 2) ref, 3) ref+signal
                'refscale', 1, ...
                'normalize', 0);
            optionNames = fieldnames(options);
            nArgs = length(varargin);
            if round(nArgs/2)~=nArgs/2
                error('propertyName/propertyValue pairs are required.')
            end
            for pair = reshape(varargin,2,[])
                if any(strncmp(pair{1}, optionNames, length(pair{1})))
                    options.(pair{1}) = pair{2};
                else
                    error('%s is not a recognized parameter name.', pair{1})
                end
            end           
                
            plotsig = bitand(options.which, 2);
            if plotsig
                if options.normalize
                    plot(this.t/options.xscale, this.yref/max(abs(this.yref)));
                else
                    plot(this.t/options.xscale, this.yref/(options.refscale*options.yscale));
                end
                hold all
            end
            
            plotsig = bitand(options.which, 1);
            if plotsig
                if options.normalize
                    plot(this.t/options.xscale, this.y/max(abs(this.y)));
                else
                    plot(this.t/options.xscale, this.y/options.yscale);
                end
            end
            
            if options.xscale == 1e-9
                xlabel('Time (ns)')
            elseif options.xscale == 1e-6
                xlabel('Time (\mu s)')
            elseif options.xscale == 1e-3
                xlabel('Time (ms)')
            elseif options.xscale ~= 1
                xlabel(['Time (' num2str(options.xscale) 's)']);
            else
                xlabel('Time (s)')
            end
            
            if options.yscale == 1e-3
                ylabel(['milli' options.yunit])
            elseif options.yscale == 1e-6
                ylabel(['\mu ' options.yunit])
            elseif options.yscale == 1e-9
                ylabel(['n ' options.yunit])
            elseif options.yscale ~= 1
                ylabel([options.yunit ' x ' num2str(options.yscale)])
            else 
                ylabel(options.yunit)
            end            
        end
        
                
        %__________________________________________________________________%
        %% ---- setting and getting properties
        function yenv = get.yenv(this)
            yenv = abs(hilbert(this.y));
        end
        
        function ydb = get.ydb(this)
            ydb = 20*log10(this.y);
        end
        
        function ny = get.ny(this)
            ny = this.y/max(abs(this.y));
        end        
        
        function yenvdb = get.yenvdb(this)
            yenvdb = 20*log10(this.yenv);
        end
        
        function nyenv = get.nyenv(this)
            nyenv = abs(hilbert(this.y/max(abs(this.y))));
        end
        
        function nydb = get.nydb(this)
            nydb = 20*log10(this.y/max(abs(this.y)));
        end
        
        function nyenvdb = get.nyenvdb(this)
            nyenvdb = 20*log10(this.nyenv);
        end        
        
        function Ts = get.Ts(this)
            Ts = mean(diff(this.t));
        end                
    end    
    
    methods(Access=protected)
        %__________________________________________________________________%
        function setOptions(this, args)
            optionNames = fieldnames(this.options);
            nArgs = length(args);
            if round(nArgs/2)~=nArgs/2
                error('propertyName/propertyValue pairs are required.')
            end

            for pair = reshape(args,2,[]) % pair is {propName;propValue}
                if any(strncmp(pair{1}, optionNames, length(pair{1})))
                    this.options.(pair{1}) = pair{2};
                else
                    error('%s is not a recognized parameter name.', pair{1})
                end
            end
        end
        
        %__________________________________________________________________%
        function isok = checkoptions(~, options, args)
            isok = 1;
            optionNames = fieldnames(options);
            nArgs = length(args);

            if round(nArgs/2)~=nArgs/2
                error('propertyName/propertyValue pairs are required.')
            end
            for pair = reshape(args,2,[])
                if any(strncmp(pair{1}, optionNames, length(pair{1})))
                    options.(pair{1}) = pair{2};
                else
                    error('%s is not a recognized parameter name.', pair{1})
                end
            end                    
        end
    end    
end

