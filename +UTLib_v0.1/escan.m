classdef escan < UTLib.utmeasurement
    %CSCAN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        E = [];
        X = [];
        step_size;
    end
    
    properties(Dependent=true,SetAccess=private)
        Edb;
    end
    
    methods
        %__________________________________________________________________%
        function this = escan(varargin)
            this@UTLib.utmeasurement(varargin{:});
        end
        
        %__________________________________________________________________%        
        function setxbystep(this, step_size, step_num)
            this.step_size = step_size;
            this.X = (0:step_num(1)-1)*step_size(1);
            this.X = this.X(:);            
        end
              
        %__________________________________________________________________%        
        % - C is a c scan 2-D array
        function setbycscan(this, C)
            [mx, r] = max(C);
            [~, c] = max(mx);
            this.E = C(r(c),:);
            this.E = this.E(:);
        end        
        
        %__________________________________________________________________%        
        % - A is an array of signals
        function setbyascanarray(this, A)
            cs = UTLib.cscan;
            cs.setcbyascanarray(A);
            this.setbycscan(cs.C);
        end               
                                
        %__________________________________________________________________%
        function normalize(this)
%             E2 = trapz(this.X,this.E.^2);
            if ~isempty(this.E)
                this.E = this.E/max(max(this.E));
%                 this.E = this.E.^2/E2;
            end
        end
        
        %__________________________________________________________________%
        function center(this)
            cent = sum(this.E.*this.X)/sum(this.E);
            [~, ind] = min(abs(this.X - cent));
            nshift = floor(length(this.X)/2 -ind);
            this.E = circshift(this.E, nshift);
        end
                
        %__________________________________________________________________%
        function plot(this, varargin)
            if ~isempty(varargin) && strcmpi(varargin{1}, 'hold')
                hold on;
            else                
                figure('units','normalized','outerposition',[0 0 1 1])
            end
            plot(this.X, this.E);
        end
    
        %__________________________________________________________________%
        function Edb = get.Edb(this)
            Edb = 20*log10(abs(this.E));
        end
        
        %__________________________________________________________________%        
        function set.E(this, e)
            this.E = e(:);
        end
        
        %__________________________________________________________________%        
        function set.X(this, x)
            this.X = abs(x(:));
        end
        
        %__________________________________________________________________%        
        function step_size = get.step_size(this)
            if isempty(this.step_size)
                step_size = mean(diff(this.X));
            else
                step_size = this.step_size;
            end
        end           
     
    end
    
end

