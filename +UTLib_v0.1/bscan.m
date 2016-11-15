classdef bscan < UTLib.utmeasurement
    %BSCAN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(SetAccess=private)
        B = [];
        X = [];
        Y = [];
    end
    
    properties(Dependent=true,SetAccess=private)
        Bdbn;
        Bn; 
        Bdb;
        Benv;
        Bnenv;
        Bnenvdb;
        Benvdb;
        ts = [];
        scan_step = [];
    end
    
    methods
        %__________________________________________________________________%
        function this = bscan(varargin)            
            this@UTLib.utmeasurement(varargin{:});            
        end    
        
        
        %__________________________________________________________________%
        function setData(this, X, Y, C, varargin)
            options = struct('tshift', 0);
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
            
            n = size(C);
            if isscalar(X) && isscalar(Y) && ~isscalar(C)                
                [this.X, this.Y] = meshgrid((0:n(2)-1)*X, (0:n(3)-1)*Y);
                if ndims(C) == 3
                    [Ct, I] = max(C,[], 3);
                    [mxc, c] = max(Ct,[],2);
                    [~, r] = max(mxc);
                    this.B = squeeze(C(r,:,:)).';
                else
                    this.B = C;
                end                 
            elseif  isvector(X) && isvector(Y) && isvector(C)
                N = [length(X), unique(Y)];
                this.X = reshape(abs(X), N)';
                this.Y = reshape(abs(Y), N)';
                this.C = reshape(C, N)';
            elseif ismatrix(X) && ismatrix(Y) && ismatrix(C)
                this.X = X;
                this.Y = Y;
                this.B = C;
                if ~isequal(size(X), size(Y), size(C))
                    error('Matrices should be equal sizes.');
                end
            else
                error('The input argument sizes or format are not correct.');
            end
            this.Y = options.tshift + this.Y;
        end
        
        
        %__________________________________________________________________%        
        % NEEDS FIXING
        function D = defectdepth(this)
            b = this.Bdbn;
            b(b>=-6) = 1;
            b(b~=1) = 0;
            meas = regionprops(b, 'BoundingBox');            
            D = [];
        end
        
        
        %__________________________________________________________________%
        function correct(this)
            Y = (this.velocity/2)*(this.Y-this.wedge_delay()*2);
            this.X = Y*sind(this.beam_angle) + this.X;
            this.Y = Y*cosd(this.beam_angle);            
        end  
        
        
        %__________________________________________________________________%
        function ind = max(this)
            [~, ind]  = max(max(this.Benv));
        end 
        
        
        %__________________________________________________________________%
        function plot(this, varargin)
            options = struct('Normalized', 1, ...
                             'Units', 'db', ...
                             'Limits', [-60 0], ....
                             'Ticks', [5 5]);
            optionNames = fieldnames(options);
            nArgs = length(varargin);
            if round(nArgs/2)~=nArgs/2
                error('propertyName/propertyValue pairs are required.')
            end

            for pair = reshape(varargin,2,[]) % pair is {propName;propValue}
                if any(strncmp(pair{1}, optionNames, length(pair{1})))
                    options.(pair{1}) = pair{2};
                else
                    error('%s is not a recognized parameter name.', pair{1})
                end
            end            
            
            if options.Normalized == 1
                if strcmpi(options.Units, 'db')
                    B = this.Bnenvdb;
                elseif strcmpi(options.Units, 'absolute')
                    B = this.Bnenv;
                end
            elseif options.Normalized == 0
                if strcmpi(options.Units, 'db')
                    B = this.Benvdb;
                elseif strcmpi(options.Units, 'absolute')
                    B = this.Benv;                    
                end                
            end
                                    
            figure('units','normalized','outerposition',[0 0 1 1],'Renderer', 'zbuffer')
            surf(this.X*1e3, this.Y*1e3, B); 
            view(2); 
            shading interp; 
            axis equal
            axis tight
            set(gca, 'FontWeight', 'normal', 'FontSize', 36, ...
                   'XTick', 0:options.Ticks(1):max(this.X(:))*1e3, ...
                   'YTick', 0:options.Ticks(2):max(this.Y(:))*1e3, ...
                   'YDir', 'reverse');
            xlabel('X (mm)')
            ylabel('Z (mm)')
            caxis(options.Limits);
%             ylim([this.thickness/2 this.thickness*3/2])
%             xlim([min(this.X(:)) max(this.X(:))]*1e3)
            
            h = colorbar;
            if strcmpi(options.Units, 'db')
                xlabel(h, 'dB', 'Units', 'normalized', 'Position', [2, -0.04]);
            elseif strcmpi(options.Units, 'absolute')
                xlabel(h, 'Units', 'Units', 'normalized', 'Position', [0.0, -0.01]);
            end   
            grid off
        end

        
        %__________________________________________________________________%
        function Bdbn = get.Bdbn(this)
            Bdbn = 20*log10(abs(this.B)/max(abs(this.B(:))));
        end    
        
        %__________________________________________________________________%
        function Bn = get.Bn(this)
            Bn = this.B/max(abs(this.B(:)));
        end            
        
        %__________________________________________________________________%
        function Bdb = get.Bdb(this)
            Bdb = 20*log10(abs(this.B));
        end     
        
        %__________________________________________________________________%
        function Benv = get.Benv(this)
            Benv = abs(hilbert(this.B));
        end     
        
        %__________________________________________________________________%
        function Bnenv = get.Bnenv(this)
            Bnenv = this.Benv/max(this.Benv(:));
        end       
        
        %__________________________________________________________________%
        function Bnenvdb = get.Bnenvdb(this)
            Bnenvdb = 20*log10(this.Bnenv);
        end               
        
        %__________________________________________________________________%
        function Bnenvdb = get.Benvdb(this)
            Bnenvdb = 20*log10(abs(this.Benv));
        end                       
        
        %__________________________________________________________________%
        function ts = get.ts(this)
           ts = this.Y(2,1) - this.Y(1,1);
        end            
        
        %__________________________________________________________________%
        function scan_step = get.scan_step(this)
            scan_step = this.X(1,2) - this.X(1,2);
        end            
    end
    
end

