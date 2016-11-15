classdef PZT < handle    
    properties (SetAccess = public, GetAccess = public)
        length = 7e-3;            % meters
        thickness = 0.2e-3;         % meters
        width = 8e-3;             % meters
        density = 7900;           % kg/m^3
        
        %coupling coefficient matrix (m/v)
        d = [0,0,-190; 0,0,-190; 0,0,450; 0,0,0; 0,0,0; 0,0,0].'*1e-12; 
        %Compliance matrix   (1/Pa)
        S = [ 1/7.6e10, 0, 0, 0, 0, 0; ... 
              0,  1/7.6e10, 0, 0, 0, 0; ...
              0, 0, 1/5.6e10, 0, 0, 0; ...
              0, 0, 0, 0, 0, 0; ...
              0, 0, 0, 0, 0, 0; ...
              0, 0, 0, 0, 0, 0];
        %permittivity matrix
        permittivity = 1850*8.854e-12*eye(3);
    end
    
    methods 
        %% Constructor
        function this = PZT(varargin)      
           
            if(size(varargin, 2) == 7)
                this.length = varargin{1};
                this.thickness = varargin{2};
                this.width = varargin{3};
                this.density = varargin{4};
                this.d = varargin{5};
                this.S = varargin{6};
                this.permittivity = varargin{7};
            elseif ~isempty(varargin{:})
                error('Wrong number of arguments entered.');
            end
        end      
        
        %% Property: set length
        function set.length(this, a)
            if ~isscalar(a) || ~isnumeric(a)
                error('The length has to be a scalar numeric value.');
            elseif a <= 0
                error('The value for the length should be positive.');
            else
                this.length = a;
            end
        end
        
        %% Property: set width 
        function set.width(this, w)
            if ~isscalar(w) || ~isnumeric(w)
                error('The width has to be a scalar numeric value.');
            elseif w <= 0
                error('The value for the width should be positive.');
            else
                this.width = w;
            end
        end
        
        %% Property: set thickness
        function set.thickness(this, tp)
            if ~isscalar(tp) || ~isnumeric(tp)
                error('The thickness has to be a scalar numeric value.');
            elseif tp <= 0
                error('The value for the thickness should be positive.');
            else
                this.thickness = tp;
            end
        end
        
        %% Property: set compliance matrix
        function set.S(this, S)            
            if ~isnumeric(S) || ~isequal(size(S),[6 6])
                error('The Young''s modulus has to be a 6x6 numeric matrix.');
            else
                this.S = S;
            end
        end
        
        %% Property: set coupling coefficient d31
        function set.d(this, d)
            if ~isnumeric(d) || ~isequal(size(d), [3 6])
                error('The coupling coefficient has to be a 3x6 numeric matrix.');
            else
                this.d = d;
            end
        end
        
        %% Property: set permittivity
        function set.permittivity(this, permittivity)
            if ~isnumeric(permittivity) || ~isequal(size(permittivity), [3 3])
                error('The permittivity should be a 3x3 numeric matrix.');
            else
                this.permittivity = permittivity;
            end
        end
        
        %% Property: set pzt density
        function set.density(this, density)
            if ~isscalar(density) || ~isnumeric(density)
                error('The PZT density has to be a scalar numeric value.');
            elseif density <= 0
                error('The value for the PZT density should be positive.');
            else
                this.density = density;
            end
        end                 
        
    end
    
end

