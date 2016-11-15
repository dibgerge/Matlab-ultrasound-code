classdef Signal < handle
    %UTSIGNAL provides operations needed for computing parameters from ultrasonic signals
    %Properties: 
    %  t: the time vector. Time should be in seconds. 
    %  y: the measured ultrasound signal. Units are arbitrary.
    %Methods: 
    %  resample:
    %  filter:
    %  removemean:
    %  gate:
    %  fft:
    %  truncate:
    %  segment:
    %  getttof:
    %  getpeaks:
    %  bandwidth:
    %  centerfrequency:
    %  plot:
    
    properties(Access=public)
        s = [];         % the main UT time domin signal
        x = [];         % the time vector corresponding to the signal        
    end
    
    properties(Dependent=true, SetAccess=private, GetAccess=public)
        Ts;             % sampling time
        Fs;             % sampling rate
        isvalid;
    end   
    
    methods
        %__________________________________________________________________%
        function self = Signal(x, s)
            if ((length(x) ~= length(s)) && ~isscalar(x))
                error(['''x'' and ''s'' should have the same length, or ''x'' ' ...
                'should be a scalar representing the sampling interval.']);            
            end
            self.s = s(:);            
            if isscalar(x)
                self.x = (0:length(s)-1)*x;
            else
                self.x = x(:);
            end
        end
        
        
        %__________________________________________________________________%        
        function s = get(self, option)
            % returns the signal according to a given option 
            if ~exist('option', 'var')
                option = '';
            end
            s = self.s;
            if ~isempty(strfind(option, 'n'))
                s = s/max(abs(s));
            end            
            if ~isempty(strfind(option, 'e'))
                s = abs(hilbert(s));
            end
            if ~isempty(strfind(option, 'd'))
                s = 20*log10(abs(s));
            end            
        end

        
        %__________________________________________________________________%
        function resample(self, ts)
            x2 = self.x(1):ts:self.x(end);            
            self.s = interp1(self.x, self.s, x2, 'spline');
            self.x = x2;
        end
        
        
        %__________________________________________________________________%
        function [f, Y] = fft(self, npts, type)
            if ~exist('npts', 'var')
                npts = [];
            end            
            if ~exist('type', 'var')
                type = '';
            end
            
            [f, Y] = signal.findfft(self.s, self.Fs, npts, type);
        end
        
        
        %__________________________________________________________________%
        function window(self, x1, x2, isindex, windfnc, varargin)
            % Apply a window to the signal within given limits
            % Inputs: 
            %   windfnc : A handle to a window function available from matlab.
            if ~exist('isindex', 'var') || isempty(isindex)
                isindex = false;
            end
            
            if ~exist('windfnc', 'var')
                windfnc = @hann;
            end
            if ~isindex
                if isempty(x1)
                    ind1 = 1;
                else
                    ind1 = find(self.x >= x1, 1);
                end
                if isempty(x2)
                    ind2 = length(self.x);
                else
                    ind2 = find(self.x >= x2, 1);
                end
            else
                if isempty(x1)
                    ind1 = 1;
                else
                    ind1 = x1;
                end
                if isempty(x2)
                    ind2 = length(self.x);
                else
                    ind2 = x2;
                end                
            end            
            if isempty(ind1) || isempty(ind2)
                error('Given time is longer than signal.');
            end
            
            N1 = ind1-1;
            N2 = ind2-ind1+1;
            N3 = length(self.x)-ind2;
            wind = windfnc(N2, varargin{:});
            wind = [zeros(N1,1); wind(:); zeros(N3,1)];
            self.s = self.s.*wind;
        end        
        
        
        %__________________________________________________________________%
        % -- filter the signal using ideal filter
        function filter(self, cutoff, type)
            if ~exist('type', 'var')
                type = 'lp';
            end
            self.s = signal.filter(self.s, self.Fs, cutoff, type);
        end        
        
        
        %__________________________________________________________________%
        % -- remove the mean of the signal
        function removemean(self)
            self.s = self.s - mean(self.s);
        end
               
        
        %__________________________________________________________________%
        function truncate(self, x1, x2, isindex)
            % Truncates a signal between the given limits x1, x2
            if ~exist('isindex', 'var')
                isindex = false;
            end            
            if ~isindex
                if isempty(x1)
                    ind1 = 1;
                else
                    ind1 = find(self.x >= x1, 1);
                end
                if isempty(x2)
                    ind2 = length(self.x);
                else
                    ind2 = find(self.x >= x2, 1);
                end
            else
                if isempty(x1)
                    ind1 = 1;
                else
                    ind1 = x1;
                end
                if isempty(x2)
                    ind2 = length(self.x);
                else
                    ind2 = x2;
                end                
            end            
            if ~isempty(ind1) && ~isempty(ind2)
                self.x = self.x(ind1:ind2);
                self.s = self.s(ind1:ind2);
            else
                error('Given time is longer than signal.');
            end
        end
        
        
        %__________________________________________________________________%
        % -- get the time of flight the ultrasound signal ---
        function varargout = gettof(this,  varargin)
            options = struct('threshold', mean(this.nyenvdb) + 13, ...
                         'holdoff', 0, ...
                         'sensitivity', 6, ...
                         'width', 2e-6, ...
                         'method', 'correlation');
            options = this.setoptions(options, varargin);
            
            us = this.segment(varargin{:});
            tof = zeros(us.N-1,1);
            if strcmpi(options.method, 'amplitude')
                I = zeros(us.N-1, 1);
                for i=1:us.N
                    [~, I(i)] = max(abs(us.ut{i}.y));
                end
                tof = diff(this.t(I));                    
            elseif strcmpi(options.method, 'envelop')
                I = zeros(us.N-1, 1);
                for i=1:us.N
                    [~, I(i)] = max(abs(hilbert(us.ut{i}.y)));
                end                
                tof = diff(this.t(I));
            elseif strcmpi(options.method, 'correlation')                
                for i=2:us.N
                    packet1 = us.ut{i-1}.y/max(abs(us.ut{i-1}.y));
                    packet2 = us.ut{i}.y/max(abs(us.ut{i}.y));
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
                varargout{2} = us;
            end
        end
        
        
        %__________________________________________________________________%
        function u = cwsegment(this, varargin)
            options = struct('start', 0, ...
                            'window', 0, ...
                            'nbins', 1);
            options = this.setoptions(options, varargin);
            startind = find(this.t >= options.start, 1);
            windind = ceil(options.window/this.Ts);
            y = repmat(this.y, 1, options.nbins);
            for i=1:options.nbins
                endind = startind+windind;
                if endind <= length(this.y)
                    y(1:startind-1,i) = 0;
                    y(endind+1:end,i) = 0;
                    y(startind:endind,i) = y(startind:endind,i).*hann(endind-startind+1);
                end
                startind = endind+1;
            end
            u = UTlib.utcollection(this.t, y);
        end
                
        
        %__________________________________________________________________%
        % -- cut the UT signal into its component reflections ---
        function seg = segment(self, varargin)
            windargs = {};
            options = struct('threshold', -10, ...
                         'holdoff', 0, ...
                         'sensitivity', 6, ...
                         'width', 2e-6, ...
                         'windfnc', @hann, ...
                         'manual' , []);
            options.windargs = {};
            options = self.setoptions(options, varargin);
            if isempty(options.manual)
                snorm = self.get('nde');
                pkloc = signal.peakfinder(snorm, options.sensitivity, options.threshold);

                % remove peaks before the holdoff time            
                pktime = self.x(pkloc);
                ind = find(pktime < options.holdoff);
                pkloc(ind) = [];    
                pktime(ind) = [];   

                % remove points within the same wave packet
                inprox = self.proximityCheck(pktime, options.width/2);
                for i=1:length(inprox)
                    % keep the maximum of all points in proximity                
                    [~, I] = max(snorm(pkloc(inprox{i}(1):inprox{i}(end))));
                    inprox{i}(I) = [];
                end            
                pkloc([inprox{:}]) = [];

                % remove end points if smaller than width 
                N = ceil(options.width/(2*self.Ts));            
                if pkloc(1) < N
                    pkloc(1) = [];
                end
                if pkloc(end)+N > length(self.s)
                    pkloc(end) = [];
                end            

                seg = cell(length(pkloc),1);
                for i=1:length(pkloc)
                    ind1 = pkloc(i)-N;
                    ind2 = pkloc(i)+N+1;
                    seg{i} = signal.Signal(self.x, self.s);
                    seg{i}.window(ind1, ind2, true, options.windfnc, options.windargs{:});
                end
            else
                seg = cell(length(options.manual)-1,1);
                for i=2:length(options.manual)
                    ind1 = options.manual(i-1);
                    ind2 = options.manual(i);
                    seg{i-1} = signal.Signal(self.x, self.s);
                    seg{i-1}.window(ind1, ind2, false, options.windfnc, options.windargs{:});
                end
            end
        end
             
        
        %__________________________________________________________________%
        function [pkind, pkmag] = getpeaks(this, varargin)            
            options = struct('sel', (max(this.yenv)-min(this.yenv))/4, ...
                'thresh', 0, ...
                'extrema', 1, ...
                'include_endpoints', 1);
            options = this.setoptions(options, varargin);
            [pkind, pkmag] = SPlib.peakfinder(this.yenv, options.sel, options.thresh, ...
                options.extrema, options.include_endpoints);
        end
        
        %__________________________________________________________________%
        function y = getfirstarrival(this, varargin)
            options = struct('threshold', -30, ...
                             'window', 1e-6, ...
                             'shift', 1e-6);
            options = this.setoptions(options, varargin);
            [ind1, ~] = this.limits(options.threshold);
            nshift = ceil(options.shift/this.Ts);
            nwind = (options.window/this.Ts);
            ind1 = ind1 - nshift;
            ind2 = ind1 + nshift +nwind-1;
            y = this.y;
%             y(1:ind1-1) = 0;
            y(ind2+1:end) = 0;
             y(ind1:ind2) = hann(nwind+nshift).*y(ind1:ind2);
        end        
        
        
        %__________________________________________________________________%        
        function [ind1, ind2] = limits(this, threshold)
            ind1 = find(this.nyenvdb > threshold, 1);
            ind2 = find(this.nyenvdb > threshold, 1, 'last');
        end

        
        %__________________________________________________________________%        
        function [bw, pbw] = bandwidth(self, threshold)
            if threshold > 0 
                error('threshold should be a dB value less than 0 dB.');
            end
            
            [f, Y] = self.fft([], 'single');
            Yn = 20*log10(abs(Y)/max(abs(Y)));
            ind = find(Yn > threshold);

            fmin = f(ind(1));
            fmax = f(ind(end));
            bw = fmax - fmin;
            pbw = 100*bw/mean([fmax, fmin]);            
        end
        
        
        %__________________________________________________________________%        
        function fc = centerfrequency(this, varargin)            
            options = struct('threshold', -12, ...
                             'type', 'mean');
            options = this.setoptions(options, varargin);
            
            if options.threshold > 0
                error('threshold should be a dB value less than 0 dB.');
            end
            
            [f, Y] = this.fft([], 'single');
            Yn = 20*log10(abs(Y)/max(abs(Y)));
            
            if strcmpi(options.type, 'mean')
                ind = find(Yn > options.threshold);    
                fmin = f(ind(1));
                fmax = f(ind(end));                
                fc = mean([fmax, fmin]);
            elseif strcmpi(options.type, 'max')
                [~,I] = max(Yn);
                fc = f(I);
            else
                error('type is not recognized. Should be ''mean'' or ''max''.');
            end
        end        
        
        
        %__________________________________________________________________%
        function r = plus(A, B)                                              
             % truncate the start
%             if A.x(1) > B.x(1)
%                 ind = find(B.x > A.x(1), 1);
%                 B.truncate(ind, [], true);
%             elseif A.x(1) < B.x(1)
%                 ind = find(A.x > B.x(1), 1);
%                 A.truncate(ind, [], true);
%             end
%             
%             % truncate the end
%             if A.x(end) > B.x(end)
%                 ind = find(A.x > B.x(end), 1);
%                 A.truncate([], ind, true);
%             elseif A.x(end) < B.x(end)
%                 ind = find(B.x > A.x(end), 1);
%                 B.truncate([], ind, true);
%             end
            
            % align the two signals            
            if A.Ts > B.Ts          
                x = A.x;
            else
                x = B.x;
            end
            y1 = interp1(A.x, A.s, x, 'spline');
            y2 = interp1(B.x, B.s, x, 'spline');
            r = signal.Signal(x, y1 + y2);
        end
        
        %__________________________________________________________________%
        function r = minus(A, B)                                              
            % truncate the start
            if A.x(1) > B.x(1)
                ind = find(B.x > A.x(1), 1);
                B.truncate(ind, [], true);
            elseif A.x(1) < B.x(1)
                ind = find(A.x > B.x(1), 1);
                A.truncate(ind, [], true);
            end
            
            % truncate the end
            if A.x(end) > B.x(end)
                ind = find(A.x > B.x(end), 1);
                A.truncate([], ind, true);
            elseif A.x(end) < B.x(end)
                ind = find(B.x > A.x(end), 1);
                B.truncate([], ind, true);
            end
            
            % align the two signals            
            if A.Ts > B.Ts          
                x = A.x;
            else
                x = B.x;
            end
            y1 = interp1(A.x, A.s, x, 'spline');
            y2 = interp1(B.x, B.s, x, 'spline');
            r = signal.Signal(x, y1 - y2);
        end
        
        
        %__________________________________________________________________%  
        function plot(self, varargin)            
            options = struct('xscale', 1e-6, ...
                'yscale', 1e-3, ...
                'yunit', 'Volts', ...
                'which', 1, ...     % two bit flag: 1) signal, 2) ref, 3) ref+signal
                'refscale', 1, ...
                'normalize', 0);            
            options = self.setoptions(options, varargin);
                
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
                    plot(self.x/options.xscale, self.s/max(abs(self.s)));
                else
                    plot(self.x/options.xscale, self.s/options.yscale);
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
        function valid = get.isvalid(self)
            valid = true;
            if ~isempty(find(isnan(self.s),1)) || ~isempty(find(isnan(self.x),1))
                valid = false;
            end                
        end
        
                
        %__________________________________________________________________%                
        function Ts = get.Ts(this)
            Ts = mean(diff(this.x));
        end 
        
        function Fs = get.Fs(this)
            Fs = 1/this.Ts;
        end         
    end    
    
    
    
    methods(Access=protected)      
        %__________________________________________________________________%
        function options = setoptions(~, options, args)
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
        
        %__________________________________________________________________%        
        function inprox = proximityCheck(~, pts, dist)
            %assume 'pts' is monotonic vector
            inprox = {};
            i = 1;
            n = 1;
            while i < length(pts)
                d = abs(pts(i:end)  - pts(i));
                ind = find(d < dist);
                ind = reshape(ind, 1, []);

                if length(ind) > 1
                    inprox{n} = i-1+ind; 
                    i = ind(end)+i;    
                    n = n+1;
                else
                    i = i + 1;
                end        
            end
        end
        
    end    
end

