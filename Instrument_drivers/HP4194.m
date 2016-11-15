classdef HP4194 < handle
    %HP4194 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        func = 1;
        hand = [];
        impedance = 2;
        gainPhase = 1;
    end
    
    methods
        % constructor
        function this = HP4194()
        end
        
        %___________________________________________________________%
        % set the function type 
        % 1: Impedance 
        % 2: gain-phase
        % 3: Impedance with probe
        function set.func(this, num)
            if isnumeric(num) && num > 0 && num < 4
                this.func = num;
                if ~isempty(this.hand) 
                    fprintf(this.hand, ['FNC' num2str(num)]);
                end
            end
        end  
        
        % Return the value of func with the corresponding description
        function val = get.func(this) 
            val.ID = this.func;
            switch(this.func)
                case 1
                    val.Name = 'Impedance';
                case 2 
                    val.Name = 'Gain-Phase';
                case 3
                    val.Name = 'Impedance with probe';
                otherwise
                    error('Invalid value found in func');
            end
        end       
        
        %___________________________________________________________%
        % set the Impedance type 
        function set.impedance(this, num)
            if isnumeric(num) && num > 0 && num < 21
                this.impedance = num;
                if ~isempty(this.hand) 
                    fprintf(this.hand, ['IMP' num2str(num)]);
                end
            end
        end
        
        % return the value of impedance with descritpion
        function val = get.impedance(this)
            val.ID = this.impedance;
            switch(this.impedance) 
                case 1
                    val.Name = '|Z|-Theta';
                case 2
                    val.Name = 'R-X';
                case 3
                    val.Name = 'Ls-Rs';
                case 4
                    val.Name = 'Ls-Q';                    
                case 5
                    val.Name = 'Cs-Rs';       
                case 6
                    val.Name = 'Cs-Q';
                case 7
                    val.Name = 'Cs-D';
                case 8
                    val.Name = '|Y|-Theta';
                case 9
                    val.Name = 'G-B';                    
                case 10
                    val.Name = 'Lp-G';                      
                case 11
                    val.Name = 'Lp-Q';
                case 12
                    val.Name = 'Cp-G';                    
                case 13
                    val.Name = 'Cp-Q';                      
                case 14
                    val.Name = 'Cp-D';
                case 15
                    val.Name = '|Z|-Ls';                    
                case 16
                    val.Name = '|Z|-Cs';                      
                case 17
                    val.Name = '|Z|-Lp';
                case 18
                    val.Name = '|Z|-Cp';
                case 19
                    val.Name = 'Lp-Rp';
                case 20
                    val.Name = 'Cp-Rp';                    
            end
        end
        
        %___________________________________________________________%
        % set the Gain-Phase type 
        function set.gainPhase(this, num)
            if isnumeric(num) && num > 0 && num < 7
                this.gainPhase = num;
                if ~isempty(this.hand) 
                    fprintf(this.hand, ['GPP' num2str(num) ';']);
                end
            end            
        end

        % return the value of gain-phase with descritpion
        function val = get.gainPhase(this)
            val.ID = this.gainPhase;
            switch(this.gainPhase) 
                case 1
                    val.Name = 'Tch/Rch [dB] - Theta';
                case 2
                    val.Name = 'Tch/Rch - Theta';
                case 3
                    val.Name = 'Tch/Rch [dB] - tau';
                case 4
                    val.Name = 'Tch/Rch [V]';                    
                case 5
                    val.Name = 'Tch/Rch [dBm]';       
                case 6
                    val.Name = 'Tch/Rch [dBV]';          
            end
        end
        
        %___________________________________________________________%        
        % set the sweeping
        function sweep(this, type, span)            
            if isnumeric(type) && type > 0  && type < 6
                cmd = ['SWP', num2str(type), ';START=', num2str(span(1)), ';STOP=',...
                    num2str(span(2)), ';NOP=' num2str(span(3))];
                 fprintf(this.hand, cmd);
            else 
                error('Connection with Impedance Analyzer not established.');
            end
        end
        
        %___________________________________________________________%        
        % set the number of averages to take
        function averaging(this, num)
            allowedVals = [1,2,4,8,16,32,64,128,256];
            if isnumeric(num) && ~isempty(find(num==allowedVals,1))
                cmd = ['NOA=' num2str(num)];
                 fprintf(this.hand, cmd);
            else 
                error('Wrong number of averages given. Should be: 1, 2, 4, 8, 16, 32, 64, 128, or 256');
            end
        end
        
        %___________________________________________________________%        
        % scale the A channel 
        function scaleR(this, scale) 
            if ~isempty(this.hand)
                if isnumeric(scale) && length(scale) == 2 
                  fprintf(this.hand, ['AMIN=' num2str(scale(1)), ';AMAX=', num2str(scale(2))]);  
                end                
            else
                error('Connection with Impedance Analyzer not established.');
            end
        end
        
        %___________________________________________________________%        
        % scale the B channel 
        function scaleI(this, scale) 
            if ~isempty(this.hand)
                if isnumeric(scale) && length(scale) == 2 
                  fprintf(this.hand, ['BMIN=' num2str(scale(1)), ';BMAX=', num2str(scale(2))]);  
                end                
            else
                error('Connection with Impedance Analyzer not established.');
            end
        end

        %___________________________________________________________%        
        % get the A channel data
        function dat = getR(this) 
            if ~isempty(this.hand)
                c = query(this.hand, 'A?');
                dat = sscanf(c, '%e,');
                dat = dat(:);
            else
                error('Connection with Impedance Analyzer not established.');
            end
        end
                
        %___________________________________________________________%
        % get the B channel data
        function dat = getI(this) 
            if ~isempty(this.hand)
                c = query(this.hand, 'B?');
                dat = sscanf(c, '%e,');
                dat = dat(:);
            else
                error('Connection with Impedance Analyzer not established.');
            end
        end
        
        %___________________________________________________________%
        % start the connection with the impedance analyzer
        function start(this)
            if isempty(this.hand)
                this.hand = instrfind('Type', 'gpib', 'BoardIndex', 0, 'PrimaryAddress', 17, 'Tag', '');

                % Create the GPIB object if it does not exist
                % otherwise use the object that was found.
                if isempty(this.hand)
                    this.hand = gpib('NI', 0, 17, 'InputBufferSize', 512000);
                else
                    fclose(this.hand);
                    this.hand = this.hand(1);
                end

                fopen(this.hand);
            end
        end
        
        % Disconnect from  the impedance analyzer
        function close(this)
            if ~isempty(this.hand)
                fclose(this.hand)
                delete(this.hand);
                this.hand = [];
            end
        end
    end
    
end

