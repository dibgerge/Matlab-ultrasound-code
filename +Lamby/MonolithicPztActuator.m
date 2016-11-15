classdef MonolithicPztActuator < Lamby.PZT
    properties
        bondThickness = 20e-6;      %meters
        bondShearModulus= 3e9;      %Pa
        V = 10;                  % Applied excitation voltage        
        plate = Lamby.IsotropicPlate();         
    end
    
    methods
        %% Constructor
        function this = MonolithicPztActuator(varargin)              
             this = this@Lamby.PZT(varargin);
        end
        
        %% Property: set bondThickness
        function set.bondThickness(this, thickness)
            if ~isscalar(thickness) || ~isnumeric(thickness) || thickness < 0
                error('The thickness has to be a positive scalar numeric value.');
            else
                this.bondThickness = thickness;
            end
        end
        
        %% Property: set bondShearModulus 
        function set.bondShearModulus(this, G)
            if ~isscalar(G) || ~isnumeric(G)
                error('The shear modulus has to be a scalar numeric value.');
            else
                this.bondShearModulus = G;
            end
        end
        
        %% Property: set voltage 
        function set.V(this, v)
            if ~isscalar(v) || ~isnumeric(v)
                error('The excitation voltage has to be a scalar numeric value.');
            else
                this.V = v;
            end            
        end
        
        %% Property: set plate 
        function set.plate(this, plate)
            import Lamby.*
            if ~isa(plate, 'IsotropicPlate')
                error('plate has to be of type IsotropicPlate.');
            else
                this.plate = plate;
            end
        end
        
        %% Method: get interfacial shear stress in bonding layer
        function [tau,tau0,x,gamma] = tau(this, type, n)
            if nargin == 1
                type = 'total';
                n = 100;
            elseif nargin == 2
                n = 100;
            end            
            
            a = this.length/2;
            Eisa = this.d(3,1)*this.V/this.thickness;
            Epzt = 1/this.S(1,1);
            phi = this.plate.Y*this.plate.thickness/(Epzt*this.thickness);
            
            if strcmpi(type, 'symmetric') || strcmpi(type, 's')
                alpha = 1; 
            elseif strcmpi(type, 'antisymmetric') || strcmpi(type, 'a')
                alpha = 3;
            elseif strcmpi(type, 'total')
                alpha = 4;
            else
                error('type is not recognized');
            end
%             dx = this.length/n;
%             x = (0:(n-1))*dx - this.length/2;
            x = linspace(-this.length/2,this.length/2,n);
            tau0 = phi*this.thickness*Epzt*Eisa/(alpha+phi);
            gamma = [];
            if this.bondThickness ~= 0 && this.bondShearModulus ~= 0
                gamma = sqrt(this.bondShearModulus*(alpha+phi)/...
                    (Epzt*this.thickness*this.bondThickness*phi));
                tau0 = tau0*gamma/cosh(gamma*a);
                tau = tau0*sinh(gamma*x);
            else
                tau = zeros(n, 1);
                tau(1) = -tau0;
                tau(end) = tau0;
            end           
        end        
        
        %% Method: get the frequency response 
        % Inputs:   f    - the frequency vector over which to compute the response
        %           mode - mode to computes the response, symmetric (s) or
        %                  antisymmetric (a)
        %           x    - (optional) The position on the plate surface where to compute response
        %                  (center of PZT assumed x = 0)
        %           sig  - (optional) If the response needs to be computed for an
        %                   arbitrary signal, sig is the given excitation
        %                   signal FFT
        %           e    - error tolerance when computing the phase velocities
        function [Hstrain, Hdisp] = frequencyResponse(this, f, mode, x, sig, threshold, e, n, space)
            import Lamby.*
            if nargin < 3
                error('not enough input arguments.');            
            end
            % Check if  f is valid and have correct values           
            if ~isnumeric(f) 
                error('Invalid values in the frequencies f. f should be positive.');
            end

            % give default value for x, if it is not given
            if nargin <= 3
                x = 0;
            elseif isempty(x)
                x = 0;                
            end
            
            if nargin <= 4 
                sig = [];
            end
            
            if ~isempty(sig)
                if length(f) ~= length(sig) 
                    error('f and sig should have the same length.');
                end 
            else
                if ~isempty(find(f<0,1))
                    error('f should be positive if not supplied with arbitrary signal.');
                end
            end
            
            if nargin <= 5  
                threshold = -30;
            elseif isempty(threshold)
                threshold = -30;
            end                            
            
            if nargin <= 6
                e = 0.1;
            end
            
            % Check if x is valid and has correct values                            
            if ~isscalar(x) || ~isnumeric(x) 
                error('x should be a scalar numeric.');
            end

            % Check if mode input is correct
            if ~strcmpi(mode, 'a') && ~strcmpi(mode, 's') && ...
                    ~strcmpi(mode, 'antisymmetric')  && ~strcmpi(mode, 'symmetric')
                error('Mode should be ''symmetric'',  ''antisymmetric'', or ''A'' or ''S''');
            end            
            
            % just get first letter of mode if it is a full word
            mode = lower(mode(1)); 
            
            if ~isempty(sig)
                % get the single sided fft
                fall = f;
                sigall = sig;
                indright = find(f>= 0);
                fright = f(indright);
                fleft = f(f < 0);
                sigright = sig(indright);            
                % Find signal limit according to the threshold
                [ind1, ind2] = findSignalLimit(sigright, threshold);
                if ind1 == 1 
                    ind1 = 2;       % ignore DC component
                end
                sig = sigright(ind1:ind2);
                f = fright(ind1:ind2);
            end
                        
            a = this.length/2;
            kl = 2*pi*f./this.plate.cl;
            kt = 2*pi*f./this.plate.ct;
            h = this.plate.thickness/2;
            f(f == 0) = 1e-6;
            [cp, fInd] = this.plate.phaseVelocity(f, mode, e);

            Hstrain = cell(length(cp), 1);
            Hdisp = cell(length(cp), 1);
            
            for i=1:length(cp)
                % equations from Giurgiutiu et al.
                k = 2*pi*f(fInd{i})./cp{i}.';
                                
                q = sqrt(kt(fInd{i}).^2-k.^2);
                p = sqrt(kl(fInd{i}).^2-k.^2);
                
                if mode == 'a'
                    Na = k.*q.*(k.^2+q.^2).*sin(p*h).*sin(q*h);
                    Da_diff = 8*k.*(k.^2-q.^2).*cos(q*h).*sin(p*h)...
                        + (k.^2-q.^2).^2.*k.*h.*(sin(p*h).*sin(q*h)./q ...
                        - cos(p*h).*cos(q*h)./p) ...
                        + (8*k.*p.*q - 4*k.^3.*(p.^2+q.^2)./(p.*q)).*sin(q*h).*cos(p*h) ...
                        + 4*k.^3.*p.*q.*h.*(-cos(p*h).*cos(q.*h)./q + sin(p*h).*sin(q.*h)./p);
                    % there was a negative sign between analytical and FEM only for asymmetric mode,
                    % so i though there some issue with the analytical
                    % solution maybe, and I add a negative sign here                    
                    F = -Na./Da_diff;
                else
                    % equations from Giurgiutiu et al.
                    Ns = k.*q.*(k.^2+q.^2).*cos(p*h).*cos(q*h);                
                    Ds_diff = 8*k.*(k.^2-q.^2).*cos(p*h).*sin(q*h)...
                        + (k.^2-q.^2).^2.*k.*h.*(sin(p*h).*sin(q*h)./p ...
                        - cos(p*h).*cos(q*h)./q) ...
                        + (8*k.*p.*q - 4*k.^3.*(p.^2+q.^2)./(p.*q)).*sin(p*h).*cos(q*h) ...
                        + 4*k.^3.*p.*q.*h.*(-cos(p*h).*cos(q.*h)./p + sin(p*h).*sin(q.*h)./q);
                    F = -Ns./Ds_diff;                    
                end
                
                [~,tau0,~,gamma] = this.tau('total');                
                if this.bondThickness == 0      %ideal bonding case
                    tau =  -tau0*sqrt(-1)*sin(k*a)./this.plate.mu;
                else                            %nonideal bonding
                    tau = tau0*sqrt(-1)*(k*sinh(gamma*a).*cos(k*a) - ...
                        gamma*cosh(gamma*a)*sin(k*a))./(k.^2 + gamma^2);
                    tau = tau/this.plate.mu;
                end
                Hstrain{i} = zeros(1,length(k));
                Hdisp{i} = zeros(1,length(k));
                if nargin > 7
                    for na=1:n+1
                        Hstrain{i} = Hstrain{i} - tau.*F.*exp(-sqrt(-1)*k*x)...
                            .*exp(-sqrt(-1)*na*(space+this.length).*k);
                        Hdisp{i} = Hdisp{i} - sqrt(-1)*tau.*F.*exp(-sqrt(-1)*k*x)...
                            .*exp(-sqrt(-1)*na*(space+this.length).*k)./k;
                    end
                else
                     Hstrain{i} = - tau.*F.*exp(-sqrt(-1)*k*x);
                     Hdisp{i} = - sqrt(-1)*tau.*F.*exp(-sqrt(-1)*k*x)./k;                    
                end
                
                % compute the response if a given arbitrary function is given
                if ~isempty(sig)                    
                    Hstrain{i} = Hstrain{i}.*sig(fInd{i});
                    Hdisp{i} = Hdisp{i}.*sig(fInd{i});
                    
                    Hstrain2{i} = zeros(1, length(sigright));
                    Hstrain2{i}(ind1+fInd{i}-1) =  Hstrain{i};
                    Hdisp2{i} = zeros(1, length(sigright));
                    Hdisp2{i}(ind1+fInd{i}-1) = Hdisp{i};  
                                
                    if mod(length(fall),2) == 1 % length of fft is odd
                        Hstrain{i} = [Hstrain2{i} conj(fliplr(Hstrain2{i}(2:end)))];
                        Hdisp{i} = [Hdisp2{i} conj(fliplr(Hdisp2{i}(2:end)))];
                    else % length of fft is even
                         Hstrain{i} = [Hstrain2{i} 0 conj(fliplr(Hstrain2{i}(2:end)))];
                         Hdisp{i} = [Hdisp2{i} 0 conj(fliplr(Hdisp2{i}(2:end)))];
                    end      
                else
                    Hstrain{i} = [zeros(1, fInd{i}(1)-1), Hstrain{i}];
                    Hdisp{i} = [zeros(1, fInd{i}(1)-1), Hdisp{i}];
                end
            end
        end
        
        %% Method: get the symmetric modes time signal response 
        % Inputs:   f    - the frequency vector over which to compute the response
        %                  corresponding to double sided FFT
        %           mode - mode to computes the response, symmetric (s) or
        %                  antisymmetric (a)
        %           x    - (optional) The position on the plate surface where to compute response
        %                  (center of PZT assumed x = 0)
        %           sig  - (optional) If the response needs to be computed for an
        %                   arbitrary signal, sig is the given excitation
        %                   signal double sided FFT
        %           threshold - sig bandwidth threshold to compute the
        %                       response. Default = -30 dB
        %           e    - error tolerance when computing the phase velocities
        function [u,e] = timeResponse(this, f, mode, x, sig, threshold, e, n, space)
            if nargin > 7
                [Hstrain, Hdisp] = this.frequencyResponse(f, ...
                    mode, x, sig, threshold, e, n, space);
            else
                [Hstrain, Hdisp] = this.frequencyResponse(f, ...
                    mode, x, sig, threshold, e);
            end
            
            nmodes = length(Hstrain);
            u = cell(nmodes,1);
            e = cell(nmodes,1);
            
            for i=1:nmodes
                u{i} = ifft(Hdisp{i});
                e{i} = ifft(Hstrain{i});
            end
        end
        
        function dispersionRemoval(this)
%                 [cga fIndga] = this.plate.antisymmetricGroupVelocity(f0,0.01);
%                 [C ind0] = min(abs(f(fInd{i})-f0));
%                 if ~isempty(ind0) 
%                     kalin = ka(ind0) + 2*pi*(f(fInd{i})-f0)./cga{i};
%                 end
%                 figure; plot(f(fInd{i})*1e-3, ka, f(fInd{i})*1e-3, kalin)
        end
        
%         function arrayFrequencyResponse(
    end  
end

