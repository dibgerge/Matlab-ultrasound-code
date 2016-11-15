classdef IsotropicPlate < handle
    %ISOTROPIC Summary of this class goes here
    %   Assume PZT actuator and sensor have the same properties
    
    properties
        % %% Plate properties
        density = 2793;    %density                    (kg/m^3)
        nu = 0.31;         %poisson's ratio
        Y = 71e9;          %Young's modulus            (Pa)
                
        thickness = 2e-3;  %Plate thickness            (m)
        
    end
    
    properties(Dependent) 
        cl;         %longitudinal wave speed    (m/s)
        ct;         %transverse wave speed      (m/s)     
        
        lambda;     %lame's constant            (Pa)
        mu;         %lame's constant            (Pa)
    end
        
    methods        
        %% Constructor
        function this = IsotropicPlate(Y, nu, density, h)
            if nargin == 4
                this.propertiesYoungPoisson(Y, nu);
                this.density = density;
                this.thickness = h;
            elseif nargin ~= 0
                error('Wrong number of arguments entered.');
            end
        end
        
        %% Property: set density
        function set.density(this, density)
            if ~isscalar(density) || ~isnumeric(density)
                error('density has to be a scalar numeric value.');
            elseif density <= 0
                error('Wrong value entered for the density.');             
            else
                this.density = density;
            end
        end
        
        %% Property: set nu
        function set.nu(this, nu)
            if ~isscalar(nu) || ~isnumeric(nu)
                error('Poisson ratio has to be a scalar numeric value.');
            else 
                this.nu  = nu;
            end
        end
        
        %% Property: set Y
        function set.Y(this, Y)
            if ~isscalar(Y) || ~isnumeric(Y)
                error('Young\''s modulus has to be a scalar numeric value.');
            elseif Y <= 0
                error('The value entered for Young\''s modulus should be positive.');
            else
                this.Y = Y;
            end
        end
        
        %% Property: set h
        function set.thickness(this, h)
            if ~isscalar(h) || ~isnumeric(h)
                error('Plate thickness has to be a scalar numeric value.');
            elseif h <= 0
                error('The value for the plate thickness should be positive.');
            else
                this.thickness = h;
            end
        end
        
        %% Property: get cl
        function cl = get.cl(this)
            cl = sqrt((this.lambda+2*this.mu)/this.density);
        end
        
        %% Property: get ct
        function ct = get.ct(this)
            ct = sqrt(this.mu/this.density);
        end
        
        %% Property: get mu 
        function mu = get.mu(this)
            mu = this.Y/(2*(1+this.nu));
        end
        
        %% Property: get lambda
        function lambda = get.lambda(this)
            lambda = this.Y*this.nu/((1+this.nu)*(1-2*this.nu));
        end                
        
        %% Method: setProperties
        function propertiesYoungPoisson(this, Y, nu)
            if(nargin == 3)
                this.Y = Y;
                this.nu = nu;
            else
                error('Wrong number of arguments entered.');
            end
        end
        
        %% Method: set Lame constants
        function propertiesLame(this, lambda, mu)
            if(nargin == 3)
               if ~isscalar(lambda) || ~isnumeric(lambda) || lambda <= 0
                  error('The lame''s constant lambda has to be a positive scalar numeric value.');
               end
               if ~isscalar(mu) || ~isnumeric(mu) || mu <= 0
                   error('The lame''s constant mu has to be a positive scalar numeric value.');
               end
               this.Y = (mu*(3*lambda+2*mu))/(mu+lambda);
               this.nu = lambda/(2*(lambda + mu));
            else
                error('Wrong number of arguments entered.');
            end
        end
        
        %% Method: set plate properties by wave speeds
        function propertiesclct(this, cl, ct)
            if(nargin == 3)
                if ~isscalar(ct) || ~isnumeric(ct) || ct <= 0
                    error('The transverse wave speed has to be a positive scalar numeric value.');
                end
                
                if ~isscalar(cl) || ~isnumeric(cl) || cl < ct
                    error('The longitudinal wave speed has be a positive scalar greater than ct.');
                end
                mu = this.density*(ct)^2;
                lambda = this.density*(cl^2-2*ct^2);
                this.propertiesLame(lambda, mu);
            else
                error('Wrong number of arguments entered.');
            end
        end    
        
        %% Method: get the phase velocity
        % Inputs:   f    - frequency scalar or vector
        %           mode - 'a' or 'antisymmetric' for antisymmetric modes
        %                  's' or 'symmetric' for symmetric modes
        %           e    - error tolerance, the smaller the more accurate
        %                  calculation. Default value = 0.1
        % %%
        % Outputs: cp - the phase velocities for the selected mode 
        %               it is a cell array, with each cell corresponding to
        %               one mode number 
        %          fInd - the indices from the given frequency vector for 
        %                 which the mode number exists
        function [cp, fInd] = phaseVelocity(this, f, mode, e)
            % Check if  2nd argument is valid and have correct values
            if ~isnumeric(f) || ~isempty(find(f<0,1))
                error('Invalid values in the frequencies f. f Should be positive.');
            end            
            
            f(f==0) = 1e-6;

            %Default value for the tolerance
            if nargin == 3
                e = 0.1;
            end

            % Check if 3rd argument is valid and has correct values
            if e <= 0 || ~isscalar(e) || ~isnumeric(e)
                error('Invalid value for the tolerance e.');
            end
            
            % Check if mode input is correct
            if ~strcmpi(mode, 'a') && ~strcmpi(mode, 's') && ...
                    ~strcmpi(mode, 'antisymmetric')  && ~strcmpi(mode, 'symmetric')
                error('Mode should be ''symmetric'',  ''antisymmetric'', or ''A'' or ''S''');
            end            
            
            % just get first letter of mode if it is a full word
            mode = lower(mode(1));                 

            c = 0:e:10000;
            cp{1} = [];
            fInd{1} = [];

            for i=1:length(f)            
                w = 2*pi*f(i);
                k2 = w.^2./c.^2;

                p = sqrt((w/this.cl)^2 - k2);
                q = sqrt((w/this.ct)^2 - k2);

                if mode == 'a'
                    % dispersion equation for antisymmetric modes
                    y = q.*tan(0.5*q*this.thickness) + ...
                        tan(0.5*p.*this.thickness).*(q.^2-k2).^2./(4.*k2.*p);
                else
                    % dispersion equation for symmetric modes
                    y = tan(0.5*q*this.thickness)./q + ...
                        4.*k2.*p.*tan(0.5*p.*this.thickness)./(q.^2-k2).^2;
                end

                d = diff(sign(y));
                %get only negative to positive zero crossings,
                %positive to negative zero crossings correspond to discontinuities
                cp_temp = c(d == 2) - e/2;

                for j=1:length(cp_temp)
                    if j > length(cp)
                        fInd{j} = i;
                        cp{j} = cp_temp(j);
                    else
                        fInd{j} = [fInd{j} i];
                        cp{j} = [cp{j} cp_temp(j)];
                    end 
                end
            end
            
            if length(cp{1}) ~= length(fInd{1})
                error('Could not find roots for some values of f, decrease the tolerance.');
            end
        end        
        
        %% Method: get the group velocity 
        function [cg, fInd] = groupVelocity(this, f, mode, e)

            % Check if f is valid and have correct values
            if ~isnumeric(f) || ~isempty(find(f<0,1))
            error('Invalid values in the frequencies f. f should be positive.');
            end

            %Default value for the tolerance
            if nargin == 3
                e = 0.1;
            end

            % Check if 3rd argument is valid and has correct values
            if e <= 0 || ~isscalar(e) || ~isnumeric(e)
                error('Invalid value for the tolerance e.');
            end
            
            % Check if mode input is correct
            if ~strcmpi(mode, 'a') && ~strcmpi(mode, 's') && ...
                    ~strcmpi(mode, 'antisymmetric')  && ~strcmpi(mode, 'symmetric')
                error('Mode should be ''symmetric'',  ''antisymmetric'', or ''A'' or ''S''');
            end        
            % just get first letter of mode if it is a full word
            mode = lower(mode(1));                 

            [cp, fInd] = this.phaseVelocity(f, mode, e);
    
            if isscalar(f)
                if f > 20e3
                    df = 1e3;
                else
                    df = f*0.1;
                end
                [cp2, fInd2] = this.phaseVelocity(f+df, mode, e);
            else
                df = f(2)-f(1);
                [cp2, fInd2] = this.phaseVelocity(f(end)+df, mode, e);
            end

            cg = cell(length(cp), 1);
            for i=1:length(cp)
                if isscalar(f)
                    cg{i} = cp{i}.^2./(cp{i} - f(fInd{i}).*(cp2{i}-cp{i})/df);
                else
                    cg{i} = cp{i}.^2./(cp{i} - f(fInd{i}).*diff([cp{i} cp2{i}])/df);    
                end        
            end
        end
        
        %% Method: ge the antisymmetric mode shapes 
        function [Ux, Uy, z] = modeshape(this, f, mode, npts)            
            % default value for npts 
            if nargin == 3
                npts = 100;
            end

            % Check if  f is valid and have correct values
            if ~isscalar(f) || ~isnumeric(f) || ~isempty(find(f<0,1))
                error('Invalid values in the frequencies f.');
            end
            
            % Check if mode input is correct
            if ~strcmpi(mode, 'a') && ~strcmpi(mode, 's') && ...
                    ~strcmpi(mode, 'antisymmetric')  && ~strcmpi(mode, 'symmetric')
                error('Mode should be ''symmetric'',  ''antisymmetric'', or ''A'' or ''S''');
            end        
            % just get first letter of mode if it is a full word
            mode = lower(mode(1));                 

            [cp, fInd] = this.phaseVelocity(f, mode);
            z = linspace(-this.thickness/2, this.thickness/2, npts);
            Ux = cell(length(cp), 1);
            Uy = cell(length(cp), 1);

            for j=1:length(cp)
                k = 2*pi*f./cp{j};
                q = 2*pi*f.*sqrt(1/this.cl^2 - 1./cp{j}.^2);
                s = 2*pi*f.*sqrt(1/this.ct^2 - 1./cp{j}.^2);
                if mode == 'a'
                    Ux{j} = (-abs((k^2-s^2))*cos(s*this.thickness/2)*sin(q*z)/...
                            (2*q*cos(q*this.thickness/2)) - s*sin(s*z)).';
                    Uy{j} = (-abs((k^2-s^2))*cos(s*this.thickness/2)*cos(q*z)/...
                            (2*k*cos(q*this.thickness/2)) + k*cos(s*z)).';
                else
                       Ux{j} = (cos(qs*z)*(ks^2-ss^2)*sin(ss*this.thickness/2)/...
                           (2*qs*sin(qs*this.thickness/2)) + ss*cos(ss*z)).';
                        Uy{j} = (-sin(qs*z)*(ks^2-ss^2)*sin(ss*this.thickness/2)/...
                            (2*ks*sin(qs*this.thickness/2)) + ks*sin(ss*z)).';     
                end
            end
        end
        
        %% Method: Phase and group velocity matching for second harmonic generation in
        %  ------- non linear guided waves
        % Arguments
        % ---------
        % type: string with values 'lame', 'longitudinal', 'crosspoint',
        % 'rayleigh'
        function [f, harmonic_mode] = velocityMatching(this, base_mode, harmonic)
            harmonic_mode = [];
            f = [];
            mode_order = str2double(base_mode(2:end));
            mode_type = lower(base_mode(1));
            
            if mode_type ~= 'a' && mode_type ~= 's'
                error('Uknown given mode type. Should be in the form <S,A><N> where N is a positive number.');
            end
            
            %-- get the Lame mode types for the given base_mode --
            if ~strcmpi(base_mode, 'a0')                
                if mode_type == 's'
                    N = mode_order*2+1;                    
                else
                    N = mode_order*2;
                end
                f(1) = N*this.ct/(this.thickness*sqrt(2));
                
                if nargin == 3
                    mult = harmonic*N;                    
                    if mod(mult,2) == 0
                        harmonic_mode{1} = ['A' num2str(mult/2)];
                    else
                        harmonic_mode{1} = ['S' num2str(floor(mult/2))];
                    end
                end
            end
            
            %-- get the longitudinal displacement symmetric modes -- 
            if mode_type ~= 'a' && ~strcmpi(base_mode, 's0')
                f(2) = mode_order*this.cl*this.ct/(sqrt(this.cl^2-this.ct^2)*this.thickness);                
                if nargin == 3
                    harmonic_mode{2} = ['S' num2str(harmonic*mode_order)];
                end
            end
            
            % -- intersections of symmetric and antisymmetric modes --
            zeta = this.ct/this.cl;
            na = 1;
            nb = 3;
            n = nb^2-na^2;
            %n = nb*(nb-1)-na*(na-1);
            %f(3) = this.ct*sqrt(n/(1-zeta^2))/(2*this.thickness);
        end
        
    end
end

