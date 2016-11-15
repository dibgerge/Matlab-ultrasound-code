classdef utnonlinear < UTLib.utsignal
    %UTNONLINEAR extracts nonlinearity parameter from ut data
    %   Detailed explanation goes here
    
    properties
        harmonic = 2;        
    end
    
    methods
        %__________________________________________________________________%
        function this = utnonlinear(t, y, varargin)
            this@UTLib.utsignal(t, y, varargin{:}); 
            
            % set the default threshold to be 6 dB from the mean of the signal
            thresh = mean(this.nyenvdb) + 6;
            % - TOFmethod : 'amplitude', 'envelop', or 'correlation'          
            %            used for computing the time of flight
            this.options = struct('threshold', thresh, ...
                         'holdoff', 0, ...
                         'width', 2e-6, ...
                         'TOFmethod', 'correlation', ...
                         'NLmethod', 'fftmax', ...
                         'NLecho', 1);
        end
        
        function beta =  getbeta(this, varargin)
            this.setOptions(varargin);
            yseg = this.segment(varargin{:});
            yseg = yseg - repmat(mean(yseg), size(yseg,1), 1);
            
            % find the frequncy of max amplitude and denote that the fundamental frequency
            if strcmpi(this.options.NLmethod, 'fftmax')                
                [~, yx] = SignalProc.filterSignal(this.t, yseg(:,this.options.NLecho), ...
                    [1e6, this.frequency*2*this.harmonic], [], 'bandpass');                
                [f, Y] = SignalProc.findfft(this.t, yx, [], 'single');                
                [A1, ind1] = max(abs(Y));
                A2 = abs(Y(find(f > this.harmonic*f(ind1), 1)));
            % integrate over a bandwidth using the given threshold
            elseif strcmpi(this.options.NLmethod, 'fftint')
                [~, yx] = SignalProc.filterSignal(this.t, yseg(:,this.options.NLecho), ...
                    [1e6, this.frequency], [], 'bandpass');
                [f, Y] = SignalProc.findfft(this.t, yx, [], 'single');
                Y = 2*Y/length(Y);
                Ydb = 20*log10(abs(Y)/max(abs(Y)));
                [~, ind] = max(Ydb);
                ind1 = ind - find(Ydb(ind:-1:1) < this.options.threshold, 1)+1;
                ind2 = ind + find(Ydb(ind:end) < this.options.threshold, 1)-1;
                A1 = trapz(f(ind1:ind2), abs(Y(ind1:ind2)));
                
                BW = f(ind2) - f(ind1);
                f1 = (f(ind2) + f(ind1))/2;                
                ind1 = find(f > f1*this.harmonic-BW/2, 1);
                ind2 = find(f > f1*this.harmonic+BW/2, 1);
                A2 = trapz(f(ind1:ind2), abs(Y(ind1:ind2)));
            % use quadrature to compute a certain frequency
            elseif strcmpi(this.options.NLmethod, 'quadrature')
                synth1 = sin(2*pi*this.frequency*this.t);
                synth2 = sin(2*pi*this.frequency*this.harmonic*this.t);
                
                A = synth1(:).*yseg(:,this.options.NLecho);
                [~, A] = SignalProc.filterSignal(this.t, A, this.frequency/2);
                [~, A] = SignalProc.findfft(this.t, A);
                A1 = abs(A(1));
                
                A = synth2(:).*yseg(:,this.options.NLecho);
                [~, A] = SignalProc.filterSignal(this.t, A, 2*this.frequency/2);
                [~, A] = SignalProc.findfft(this.t, A);                
                A2 = abs(A(1));
            end
            beta = A2./A1.^2;            
        end
        
    end
    
end

