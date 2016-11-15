classdef utcollection < handle
    %UTCOLLECTION performs operations on an array of utsignal objects
    %   Detailed explanation goes here
    
    properties
        ut;         % the utsignals
    end
    
    properties(Dependent=true, SetAccess=private, GetAccess=public)
        N;       % number of signals
    end
    
    methods
        %__________________________________________________________________%
        function this = utcollection(varargin)
            if mod(nargin,2) ~= 0
                error('Constructor inputs should be of the form t1, y1, t2, y2, ...');
            elseif nargin == 2 && size(varargin{2}, 2) > 1
                N = size(varargin{2},2);
                this.ut = cell(N,1);
                for i=1:N
                    this.ut{i} = UTlib.utsignal(varargin{1}, varargin{2}(:,i));
                end                
            else
                N = nargin/2;
                this.ut = cell(N,1);
                for i=1:N
                    this.ut{i} = UTlib.utsignal(varargin{2*i-1}, varargin{2*i});
                end                
            end
        end

        
        %__________________________________________________________________%
        function add(this, varargin)
            if mod(nargin-1,2) ~= 0
                error('Constructor inputs should be of the form t1, y1, t2, y2, ...');
            end
            N = (nargin-1)/2;
            for i=1:N
                this.ut{end+i} = UTlib.utsignal(varargin{2*i-1}, varargin{2*i});
            end
        end

        
        %__________________________________________________________________% 
        function removeinitbang(this, varargin)
             for i=1:this.N
                 this.ut{i}.removeinitbang(varargin{:});
             end            
        end
        
        %__________________________________________________________________%        
         function truncate(this, varargin)
             for i=1:this.N
                 this.ut{i}.truncate(varargin{:});
             end
         end
         
        %__________________________________________________________________%
        function [f, Y] = fft(this, npts, type)
            this.synchronize(); 
            if nargin == 1
                npts = [];
                type = 'shift';                
            elseif nargin == 2
                type = 'shift';
            end
            for i=1:this.N
                [f, Y(:,i)] = this.ut{i}.fft(npts, type);
            end
        end
         
        
        %__________________________________________________________________%
        function filter(this, varargin)
            for i=1:this.N
                this.ut{i}.filter(varargin{:});
            end            
        end          
        
        
        %__________________________________________________________________%
        function yseg = segment(this, varargin)
            yseg = cell(this.N,1);
            for i=1:this.N
                yseg{i} = this.ut{i}.segment(varargin{:});
            end            
        end        
        

        %__________________________________________________________________%
        function synchronize(this)
            % truncate the start and end
            tstart = -Inf;
            tend = Inf;
            ts = Inf;
            for i=1:this.N
                if this.ut{i}.t(1) > tstart
                    tstart = this.ut{i}.t(1);
                end
                if this.ut{i}.t(end) < tend
                    tend = this.ut{i}.t(end);
                end
                if this.ut{i}.Ts < ts;
                    ts = this.ut{i}.Ts;
                end
            end
            
            tnew = (tstart:ts:tend)';
            for i=1:this.N
                this.ut{i}.y = interp1(this.ut{i}.t, this.ut{i}.y, tnew, 'spline');
                this.ut{i}.t = tnew;
            end                
        end
        
        
        %__________________________________________________________________%
        function plot(this)
            for i=1:this.N
                this.ut{i}.plot()
                hold all;
            end
        end               
        
        
        %__________________________________________________________________%        
        function n = get.N(this)
            n = length(this.ut);
        end        
        
    end
    
end

