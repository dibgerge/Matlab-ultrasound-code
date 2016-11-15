function [t, varargout] = syncSignals(varargin)

% --- Default options values --- %
options = struct('dt', 0, ...
                 'T', 0);

%# read the acceptable names
optionNames = fieldnames(options);

%# count arguments
nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2
   error('syncSignals needs propertyName/propertyValue pairs.')
end

n = 1;
signals = struct('t', [], 'data', []);
for pair = reshape(varargin,2,[]) %# pair is {propName;propValue}
   if ~ischar(pair{1})
       signals(n).t = pair{1};
       signals(n).data = pair{2};
       n=n+1;
   else
    if any(strmatch(pair{1}, optionNames))
          %# overwrite options. If you want you can test for the right class here
        %# Also, if you find out that there is an option you keep getting wrong,
        %# you can use "if strcmp(inpName,'problemOption'),testMore,end"-statements
        options.(pair{1}) = pair{2};
    else
      error('%s is not a recognized parameter name',pair{1})
    end
   end   
end

%%________________________________________________________________________%

if options.dt == 0    
    options.dt = Inf;
    for i=1:length(signals)
        dt = mean(diff(signals(i).t));
        if  dt < options.dt
            options.dt = dt;
        end
    end
end

if options.T == 0
    for i=1:length(signals)
        T = signals(i).t(end);
        if  T > options.T
            options.T = T*2;
        end
    end
end

t = (0:options.dt:options.T)';
varargout = cell(length(signals),1);
for i=1:length(signals)
    varargout{i} = interp1(signals(i).t, signals(i).data, t, 'linear', 0);
end

end

