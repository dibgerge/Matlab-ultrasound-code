classdef gwnonlinear < handle
    %UTNONLINEARGW extracts nonlinearity parameters fro nonlinear guided
    %wave measurements
    %   Detailed explanation goes here
    
    properties
        t = [];
        y = [];
        frequency = 1.1e6;
        harmonic = 2;
    end
    
    properties(Access=protected)
        options;
        S = [];
        F = []; 
        T = []; 
        P = [];
    end    
    
    properties(Dependent=true)
        Ts;
    end        
    
    methods
        %__________________________________________________________________%
        function this = gwnonlinear(t, y, varargin)
            this.t = t;
            this.y = y;
                        
            this.options = struct('window', 128, ...
                         'overlap', 127, 'threshold', -10, 'gate', 20e-6);
            this.setOptions(varargin);
        end
        
        %__________________________________________________________________%        
        function [S, F, T, P] = spectrogram(this, varargin)
            this.setOptions(varargin);
            [S, F, T, P] = spectrogram(this.y, this.options.window, this.options.overlap, ...
                    this.options.window*2, 1/this.Ts);
            this.S = S;
            this.F = F;
            this.T = T;
            this.P = P;
        end
        
        %__________________________________________________________________%        
        function [T, y1] = specgramatf(this, varargin)
            if ~isempty(varargin) || isempty(this.S)
                this.setOptions(varargin);            
                this.spectrogram(varargin{:});
            end
            indf1 = find(this.F>this.frequency, 1);
            y1 = this.S(indf1,:);
            T = this.T;
        end        
        
        %__________________________________________________________________%        
        function [T, y1] = specgramathr(this, varargin)
            if ~isempty(varargin) || isempty(this.S)
                this.setOptions(varargin);
                this.spectrogram(varargin{:});
            end
            indf1 = find(this.F>this.frequency*this.harmonic, 1);
            y1 = this.S(indf1,:);
            T = this.T;
        end                

        %__________________________________________________________________%        
        function b = beta(this, varargin)
            [T1, y1] = this.specgramatf(varargin{:});
            [T2, y2] = this.specgramathr();
            
            startInd = find(20*log10(y1/max(y1)) > this.options.threshold, 1);
            endInd = find(T1 > T1(startInd) + this.options.gate, 1);
            A1 = trapz(T1(startInd:endInd), abs(y1(startInd:endInd)));


            startInd = find(20*log10(y2/max(y2)) > this.options.threshold, 1);
            endInd = find(T2 > T2(startInd) + this.options.gate, 1);
            A2 = trapz(T2(startInd:endInd), abs(y2(startInd:endInd)));    
            b = A2/A1.^2;
        end                
        
        %__________________________________________________________________%                
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
    end      
    
end

