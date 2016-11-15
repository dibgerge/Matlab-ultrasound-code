classdef Material < handle
    %MATERIAL Specify the properties of homogeneous isotropic materials
    %Properties:
    %  density  : in kg/m^3
    %  cl       : in m/s
    %  ct       : in m/s
    % the properties defaults are for Plexiglass.
    %Methods:
    %  none
    
    properties
        density;    % in kg/m^3
        cl;         % in m/s
        ct;         % in m/s
    end
    
    methods
        %__________________________________________________________________%        
        function this = Material(varargin)
             % Defaults for the required inputs
            options = struct('density', 1180, ...
                             'cl', 2680, ...
                             'ct', 1320);
                         
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
            this.density = options.density;
            this.cl = options.cl;
            this.ct = options.ct;
        end
        
        %__________________________________________________________________%
        function set.density(this, val)
            if isscalar(val) && val > 0
                this.density = val;
            else
                error('density should be a positive scalar.');
            end
        end        
        
        %__________________________________________________________________%
        function set.cl(this, val)
            if isscalar(val) && val > 0
                this.cl = val;
            else
                error('cl should be a positive scalar.');
            end            
        end
        
        %__________________________________________________________________%
        function set.ct(this, val)
            if isscalar(val) && val > 0
                this.ct = val;
            else
                error('ct should be a positive scalar.');
            end
        end                      
    end
    
end

