classdef Inspection < handle
    %INSPECITON Contains description and methods needed to define an
    %  inspection
    %Properties:
    %  wedge
    
    properties
        wedge;
        transducer;
        material;
        thickness;
    end
    
    properties(Dependent=true,SetAccess=private)
        mode;
        beam_angle;
        lambda;
        near_field;
    end
    
    methods
        %__________________________________________________________________%        
        function this = Inspection(varargin)
            % Defaults for the required inputs
            options = struct('wedge', UTlib.Wedge(), ...
                             'transducer', UTlib.Transducer(), ...
                             'material', UTlib.Material(), ...
                             'thickness', 1e-3);

            % Make sure that the inputs are valid
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
            this.wedge = options.wedge;
            this.transducer = options.transducer;
            this.material = options.material;
            this.thickness = options.thickness;
        end
        
        
        %__________________________________________________________________%
        function sp = beamspread(this, limit)
            if limit == -Inf
                c = 1.22;
            elseif limit == -6
                c = 0.56;
            elseif limit == -20
                c = 1.08;
            else
                error('limit should have a numeric value of -Inf dB, -6 dB, or -20 dB.');
            end                
            sp = 2*asind(c*this.lambda/this.transducer.diameter);
        end
        
                
        %__________________________________________________________________%
        % Get functions for the dependent properties        
        %__________________________________________________________________%
        
        %__________________________________________________________________%
        function lambda = get.lambda(this)
            lambda = [];
            if ~isempty(strfind(this.mode, 'L'))
                lambda = [lambda, this.material.cl/this.transducer.centerfrequency()];
            end
            
            if ~isempty(strfind(this.mode, 'S'))
                lambda = [lambda, this.material.ct/this.transducer.centerfrequency()];
            end
        end
        
        
        %__________________________________________________________________%
        function mode = get.mode(this)
            % application of snell's law
            % assume longitudinal wave is propagating in the wedge
            rl = this.material.cl*sind(this.wedge.angle)/this.wedge.material.cl;
            rt = this.material.ct*sind(this.wedge.angle)/this.wedge.material.cl;
            
            if rl > 1 && rt > 1
                mode = '';
            elseif rl > 1 && rt <= 1
                mode = 'S';
            elseif rl <= 1 && rt > 1
                mode = 'L';
            elseif rl <= 1 && rt <= 1
                mode = 'LS';
            end
        end
        
        
        %__________________________________________________________________%
        function nf = get.near_field(this)
               nf = this.transducer.diameter^2./(4*this.lambda);
        end
        
        
        %__________________________________________________________________%
        function beamangle = get.beam_angle(this)            
            mode = this.mode;
            beamangle = [];
            if ~isempty(strfind(mode, 'L'))
                beamangle = [beamangle, ...
                    asind(this.material.cl*sind(this.wedge.angle)/this.wedge.material.cl)];
            end
            
            if ~isempty(strfind(mode, 'S'))
                beamangle = [beamangle, ...
                    asind(this.material.ct*sind(this.wedge.angle)/this.wedge.material.cl)];
            end
        end        
        
        %__________________________________________________________________%
        % Set functions for properties        
        %__________________________________________________________________%
        
        %__________________________________________________________________%
        function set.wedge(this, val)
            if isa(val, 'UTlib.Wedge')
                this.wedge = val;
            else
                error('wedge should be of type UTlib.Wedge().');
            end
        end
        
        
        %__________________________________________________________________%
        function set.transducer(this, val)
            if isa(val, 'UTlib.Transducer')
                this.transducer = val;
            else
                error('transducer should be of type UTlib.Transducer().');
            end
        end
        
        
        %__________________________________________________________________%
        function set.material(this, val)
            if isa(val, 'UTlib.Material')
                this.material = val;
            else
                error('material should be of type UTlib.Material().');
            end
        end
        
        
        %__________________________________________________________________%
        function set.thickness(this, val)
            if isscalar(val) && val >= 0
                this.thickness = val;
            else
                error('thickness should be scalar and positive.');
            end
        end
    end
    
end

