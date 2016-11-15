classdef Transducer < handle
    %TRANSDUCER defines a transducer response, and its corresponding wedge, if any
    %Properties: 
    %  response : utsignal object, storing the response of this transducer
    %             off of a side drilled hole. 
    %  diameter : in meters, the diameter of the transducer element.
    %Dependent properties:
    %  beam_angle: 
    
    properties
        response;
        diameter;
    end
    
    methods
        %__________________________________________________________________%
        function this = Transducer(varargin)
            % Defaults for the required inputs
            options = struct('response', UTlib.utsignal(0, 0), ...
                             'diameter', 6.35e-3);
                         
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
            this.response = options.response;
            this.diameter = options.diameter;
        end

        
        %__________________________________________________________________%
        function [bw, pbw] = bandwidth(this, threshold)
            [bw, pbw] = this.response.bandwidth(threshold);
        end
        
        %__________________________________________________________________%
        function fc = centerfrequency(this, type)            
            if nargin == 1
                type = 'max';
            end                
            fc = this.response.centerfrequency('type', type);
        end
    end
    
end

