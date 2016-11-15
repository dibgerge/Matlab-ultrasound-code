classdef Bscan < handle
    %BSCAN Defines an ultrasound  bscan
    %   Detailed explanation goes here
    
    properties(SetAccess=private)
        inspection;
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
        Ts = [];
        scan_step = [];
    end
    
    methods
        %__________________________________________________________________%
        function this = Bscan(inspection)            
            if nargin == 1
                if isa(inspection, 'UTlib.Inspection')
                    this.inspection = inspection;
                end
            end
        end    
        
        
        %__________________________________________________________________%
        function setdata(this, X, Y, C, varargin)
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
            
            if ~ismatrix(C) || length(n) ~= 2
                error('C should be a matrix.');
            end
            
            if isscalar(X) && isscalar(Y)
                [this.X, this.Y] = meshgrid((0:n(2)-1)*X, (0:n(1)-1)*Y);
                this.B = C;                
            elseif ismatrix(X) && ismatrix(Y)
                if ~isequal(size(X), size(Y), size(C))
                    error('Matrices should be equal sizes.');
                end
                this.X = X;
                this.Y = Y;
                this.B = C;                
            else
                error('The input argument sizes or format are not correct.');
            end
            this.Y = options.tshift + this.Y;
        end
        
        
        %__________________________________________________________________%        
        function mx = tipamp(this, thresh)            
            Bthresh = this.Bnenvdb;
            Bthresh(Bthresh > thresh) = 1;
            Bthresh(Bthresh < 1) = 0;
            cc = bwconncomp(Bthresh);
            stats = regionprops(cc, 'Centroid');
            Nr = length(cc.PixelIdxList);            
            regionmaxes = zeros(Nr, 1);           
            for i=1:Nr
                regionmaxes(i) =  max(this.Bnenvdb(cc.PixelIdxList{i}));
            end            
            [~, theregion] = max(regionmaxes);
            for i=1:Nr
                if stats(i).Centroid(2) < stats(theregion).Centroid(2)
                    theregion = i;
                end
            end                        
            mx = regionmaxes(theregion);
        end
        
        
        %__________________________________________________________________%
        function correct(this)            
            if strcmpi(this.inspection.mode, 'S')
                c = this.inspection.material.ct;
            elseif strcmpi(this.inspection.mode, 'L')
                c = this.inspection.material.ct;
            else
                error('The number of modes in the structure is more than one, or do not exist.')
            end
                            
            Y = (c/2)*(this.Y-this.inspection.wedge.delay*2);
            this.X = Y*sind(this.inspection.beam_angle) + this.X;
            this.Y = Y*cosd(this.inspection.beam_angle);
        end  
        
        
        %__________________________________________________________________%
        function [mx, I, J] = max(this)
            [mx, ind]  = max(this.B(:));
            [I, J] = ind2sub(size(this.B), ind);
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
            set(gca, 'FontWeight', 'normal', 'FontSize', 36, 'YDir', 'reverse', ...
                    'XTick', 0:options.Ticks(1):max(this.X(:))*1e3, ...
                    'YTick', 0:options.Ticks(2):max(this.Y(:))*1e3);   
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
        function Ts = get.Ts(this)
           Ts = this.Y(2,1) - this.Y(1,1);
        end            
        
        %__________________________________________________________________%
        function scan_step = get.scan_step(this)
            scan_step = this.X(1,2) - this.X(1,1);
        end            
    end
    
end

