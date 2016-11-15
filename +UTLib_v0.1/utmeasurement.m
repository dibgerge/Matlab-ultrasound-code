classdef utmeasurement < handle
    %UTMEASUREMENT represents the measurement transducer and sample used
    %   Detailed explanation goes here
    
    properties
        trans_dia = [];         % meters
        thickness = [];         % meters
        width = [];
        length = []
        frequency = [];         % Hz
        velocity = [];          % m/s
        wedge_len = [];
        wedge_speed = [];
        beam_angle = [];
        scan_dir = '+';
        mode = '';
    end
    
    properties(Dependent=true)
        lambda = [];
        k = [];
        wedge_delay = [];
    end
    
    methods
        %__________________________________________________________________%
        function this = utmeasurement(varargin)
            options = struct('trans_dia', 0.25*25.4e-3, ...
                        'thickness', 5e-3, ...
                        'width', [], ...
                        'length', [], ...
                        'frequency', 5e6, ...
                        'wedge_len', [], ...
                        'wedge_speed', [], ...
                        'mode', 'SV', ...
                        'velocity', [], ...
                        'beam_angle', 90, ...
                        'scan_dir', '+');
            % case of a copy constructor
            if length(varargin) == 1 && isa(varargin{1}, 'UTLib.utmeasurement')
                this.trans_dia = varargin{1}.trans_dia;
                this.thickness = varargin{1}.thickness;
                this.frequency = varargin{1}.frequency;
                this.velocity = varargin{1}.velocity;            
                this.wedge_len = varargin{1}.wedge_len;
                this.wedge_speed = varargin{1}.wedge_speed;
                this.beam_angle = varargin{1}.beam_angle;
                this.scan_dir = varargin{1}.scan_dir;
                this.mode = varargin{1}.mode;
            else
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
                this.trans_dia = options.trans_dia;
                this.thickness = options.thickness;
                this.frequency = options.frequency;
                this.velocity = options.velocity;            
                this.wedge_len = options.wedge_len;
                this.wedge_speed = options.wedge_speed;
                this.beam_angle = options.beam_angle;
                this.scan_dir = options.scan_dir;
                this.mode = options.mode;
            end
        end
        
        %__________________________________________________________________%
        function nf = getnearfield(this)            
            nf = this.trans_dia^2/(4*this.lambda);
        end

        %__________________________________________________________________%
        function beam = beamprofile(this, D)         
            beam = 1 - exp(sqrt(-1)*this.k*this.trans_dia.^2./D);
        end
        
        %__________________________________________________________________%
        function Dhat = diffraction(this, D)
            u = this.k*(this.trans_dia^2)./D;
            Dhat = 1 - exp(sqrt(-1)*u).*(besselj(0, u) - sqrt(-1)*besselj(1, u));
        end        
        
        %__________________________________________________________________%        
        % The -inf dB beam spread in degrees
        function sp = getspread(this)  
            sp = 2*asind(1.22*this.lambda/this.trans_dia);
        end
        
        % Dependent properties
        %__________________________________________________________________% 
        function lambda = get.lambda(this)
            lambda = this.velocity./this.frequency;
        end
        
        function k = get.k(this)
            k = 2*pi./this.lambda;
        end        
        
        function wedge_delay = get.wedge_delay(this)
            wedge_delay = this.wedge_len/this.wedge_speed;
        end
        
    end
    
end

