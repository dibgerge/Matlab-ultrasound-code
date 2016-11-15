classdef Cscan < handle
    %CSCAN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(SetAccess=private)
        inspection;
        C = [];
        X = [];
        Y = [];
    end
    
    properties(Dependent=true,SetAccess=private)
        Cdb;
        Cdbn;
        Cn;
        Xstep;
        Ystep;
    end
    
    methods
        %__________________________________________________________________%
        function this = Cscan(inspection)            
             if nargin == 1
                if isa(inspection, 'UTlib.Inspection')
                    this.inspection = inspection;
                end
            end          
        end
        
        
        %__________________________________________________________________%                
        function setdata(this, X, Y, C)
            n = size(C);
            if isscalar(X) && isscalar(Y) && ~isscalar(C)
                [this.X, this.Y] = meshgrid((0:n(2)-1)*X, (0:n(1)-1)*Y);                 
                if ndims(C) == 3
                    for i=1:n(1)
                        for j=1:n(2)
                            this.C(i,j) = max(abs(hilbert(squeeze(C(i,j,:)))));
                        end
                    end
                else
                    this.C = C;
                end                    
            elseif  isvector(X) && isvector(Y) && isvector(C)
                N = [length(unique(X)), length(unique(Y))];
                this.X = reshape(abs(X), N)';
                this.Y = reshape(abs(Y), N)';
                this.C = reshape(C, N)';
                % determine if flipping is required do to a raster scan 
                % Assuming that X is always the direction of the scan 
                if isempty(find(this.X(1,:) - fliplr(this.X(2,:)),1))
                    this.X(2:2:end,:) = fliplr(this.X(2:2:end,:)); 
                    this.C(2:2:end,:) = fliplr(this.C(2:2:end,:));
                end                                    
            elseif ismatrix(X) && ismatrix(Y) && ismatrix(C)
                this.X = X;
                this.Y = Y;
                this.C = C;
                if ~isequal(size(X), size(Y), size(C(:,:,1)))
                    error('Matrices should be equal sizes.');
                end
                if ndims(C) == 3
                    for i=1:n(1)
                        for j=1:n(2)
                            this.C(i,j) = max(abs(hilbert(squeeze(C(i,j,:)))));
                        end
                    end
                end
            else
                error('The input argument sizes or format are not correct.');
            end           
        end
                        
        
        %__________________________________________________________________%        
        function [xind, yind, val] = max(this)
            [mxy, yind] = max(this.C);
            [~, xind] = max(mxy);
            yind = yind(xind);
            val = this.C(yind, xind);
        end

        
        %__________________________________________________________________%
        function center(this, thresh)            
            if nargin == 1;
                thresh = mean(this.Cdb(:)) + 6;
            end

            binaryImage = this.Cdb > thresh;

            measurements = regionprops(binaryImage, 'Centroid');
            cent = measurements.Centroid;

            [r, c] = size(binaryImage);
            rowsToShift = round(r/2- cent(2));
            columnsToShift = round(c/2 - cent(1));

            % Call circshift to move region to the center.
            this.C = circshift(this.C, [rowsToShift columnsToShift]);
        end
        
        
        %__________________________________________________________________%
        function medfilt(this, order)
            this.C = medfilt2(this.C, order);            
        end
        
        
        %__________________________________________________________________%        
        function L = defectlength(this, dir)
            c = this.Cdbn;
            c(c>=-6) = 1;
            c(c~=1) = 0;
            meas = regionprops(c, 'BoundingBox');
            if dir == 'x'
                L = meas.BoundingBox(3)*this.step_size(1);
            elseif dir == 'y'
                L = meas.BoundingBox(4)*this.step_size(2);
            end
        end
        
        
        %__________________________________________________________________%
        function plot(this, varargin)
            options = struct('Normalized', 1, ...
                             'Units', 'db', ...
                             'Limits', [-20 0], ....
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
                    C = this.Cdbn;
                elseif strcmpi(options.Units, 'absolute')
                    C = this.Cn;
                end
            elseif options.Normalized == 0
                if strcmpi(options.Units, 'db')
                    C = this.Cdb;
                elseif strcmpi(options.Units, 'absolute')
                    C = this.C;
                    
                end                
            end
            figure('units','normalized','outerposition',[0 0 1 1],'Renderer', 'painters')
            surf(this.X*1e3, this.Y*1e3, C)
            view(2)
            shading interp
            axis equal
            axis tight
            set(gca, 'FontWeight', 'normal', 'FontSize', 36, ...
                   'XTick', 0:options.Ticks(1):max(1e3*this.X(:)), ...
                   'YTick', 0:options.Ticks(2):max(1e3*this.Y(:)));
            xlabel('X (mm)')
            ylabel('Y (mm)')
            caxis(options.Limits);
            h = colorbar;
            if strcmpi(options.Units, 'db')
                xlabel(h, 'dB', 'Units', 'normalized', 'Position', [2, -0.03]);
            elseif strcmpi(options.Units, 'absolute')
                xlabel(h, 'Units', 'Units', 'normalized', 'Position', [0.0, -0.01]);
            end
        end
    
        %__________________________________________________________________%
        function Cdb = get.Cdb(this)
            Cdb = 20*log10(abs(this.C));
        end
        
        %__________________________________________________________________%
        function Cdbn = get.Cdbn(this)
            Cdbn = 20*log10(abs(this.C)/max(abs(this.C(:))));
        end        
        
        %__________________________________________________________________%
        function Cn = get.Cn(this)
            if ~isempty(this.C)
                Cn = this.C/max(max(this.C));
            end
        end
        
    end
    
end

