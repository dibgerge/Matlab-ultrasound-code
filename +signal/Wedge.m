classdef Wedge < handle
    %WEDGE defines the dimensions and materials of wedge for angle
    %  beam inspection
    %   
    %   
    %         ___________________
    %        /                  |
    %      /\                   |
    %    /   \ L3               |
    %  /angle \                 |
    %/_________\________________|    
    %     L2           L1
    %Properties: 
    %  L1: in meters
    %  L2: in meters
    %  L3: in meters
    %  L4: in meters
    %  angle: in degrees
    %  material: a class of type Modelib.Material
    %  delay: in second is a dependent property
    %  exists: uses the specified dimensions of the wedge to determine if
    %          they are physically feasible.
    %Methods:
    %  none
    
    properties
        L1;     % in m
        L2;     % in m
        L3;     % in m
        L4;     % in m
        angle;  % in Degrees
        material;
    end
    
    properties(Dependent=true,SetAccess=private)
        delay;  % in seconds
        exist;
    end
    
    methods
        %__________________________________________________________________%
        function this = Wedge(varargin)
            % Defaults for the required inputs
            options = struct('L1', 0, ...
                             'L2', 0, ...
                             'L3', 0, ...
                             'L4', 0, ...
                             'angle', 0, ...
                             'material', UTlib.Material());
                         
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
            this.L1 = options.L1;
            this.L2 = options.L2;
            this.L3 = options.L3;
            this.L4 = options.L4;
            this.angle = options.angle;
            this.material = options.material;
        end
        
        %__________________________________________________________________%
        function set.L1(this, val)
            if isscalar(val) && val >= 0
                this.L1 = val;
            else
                error('L1 should be a positive scalar.');
            end
        end
        
        %__________________________________________________________________%
        function set.L2(this, val)
            if isscalar(val) && val >= 0
                this.L2 = val;
            else
                error('L2 should be a positive scalar.');
            end
        end

        %__________________________________________________________________%
        function set.L3(this, val)
            if isscalar(val) && val >= 0
                this.L3 = val;
            else
                error('L3 should be a positive scalar.');
            end
        end
        
        %__________________________________________________________________%
        function set.L4(this, val)
            if isscalar(val) && val >= 0
                this.L4 = val;
            else
                error('L4 should be a positive scalar.');
            end
        end

        %__________________________________________________________________%
        % the angle of the wedge with respect to the vertical.
        function set.angle(this, val)
            if isscalar(val) && val >= 0
                this.angle = val;
            else
                error('angle should be a positive scalar.');
            end
        end        

        %__________________________________________________________________%
        function set.material(this, val)
            if isa(val, 'UTlib.Material')
                this.material = val;
            else
                error('material should be of type UTlib.Material.');
            end
        end        

        %__________________________________________________________________%
        function delay = get.delay(this)
            % Assume the longitudinal wave is excited in the wedge always
            delay = this.L3/this.material.cl;            
        end
        
        %__________________________________________________________________%
        function exist = get.exist(this)
            exist = false;
            if this.L1 > 0 && this.L2 > 0 && this.L3 > 0 && this.L4 > 0
                exist = true;
            end                
        end            
                
        
    end
    
end
