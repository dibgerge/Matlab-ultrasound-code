function [pod, plow, phigh, ai, ahat, p_ahat, threshold, pf] = findPOD6(F, varargin)

% --- Default options values --- %
options = struct('PDFParams', [0 1], ...
                 'maxDefect', 10e-3, ...
                 'xiRange', [0 2], ...
                 'NoiseLevel', 0,...
                 'FArate', 0.05,...
                 'aPts', 100, ...
                 'ahatPts', 100,...
                 'Threshold', [], ...
                 'PD', [], ...
                 'samplePts', 1e6);

% --- read the acceptable names ---%
optionNames = fieldnames(options);

% --- count arguments ---%
nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2
   error('findPOD needs propertyName/propertyValue pairs.')
end

for pair = reshape(varargin,2,[]) % pair is {propName;propValue}
   if any(strmatch(pair{1}, optionNames))
      % overwrite options. If you want you can test for the right class here
      % Also, if you find out that there is an option you keep getting wrong,
      % you can use "if strcmp(inpName,'problemOption'),testMore,end"-statements
      options.(pair{1}) = pair{2};
   else
      error('%s is not a recognized parameter name',pair{1})
   end
end
 
%% --- (1) Initialization --- %
ai = linspace(0, options.maxDefect, options.aPts);
p_ahat = cell(options.aPts,1);
ahat = cell(options.aPts,1);
% p_ahat = zeros(options.ahatPts, options.aPts);
% ahat = zeros(options.ahatPts, options.aPts);

%% --- (2) Obtain the variability PDF ---%
if options.PDFParams(2) > 0
    p_xi = makedist('Normal', 'mu', options.PDFParams(1), 'sigma', options.PDFParams(2));
    p_xi = p_xi.truncate(options.xiRange(1), options.xiRange(end));
end

for i=1:options.aPts    
    %% --- (3) Get samples and find the 'ahat' PDF  ---%
    if options.PDFParams(2) > 0
        xivals = random(p_xi, options.samplePts, 1);
    
        [xi, a] = ndgrid(xivals, ai(i));
        y = F(xi, a);
    else
        xi = ones(options.samplePts,1)*options.PDFParams(1);
        a = ai(i)*ones(options.samplePts,1);
        y = F(xi, a);        
    end          
    
    %% --- (4) Addd noise to the problem --- %
    if options.NoiseLevel > 0        
        noisePDF = makedist('Normal', 'mu', 0, 'sigma', options.NoiseLevel);
%         noisePDF = noisePDF.truncate(0, inf);

        noise = random(noisePDF, options.samplePts, 1);
        y = abs(y + abs(hilbert(noise)));
    end
    
    [p, val] = hist(y, options.ahatPts); 
    ahat{i}= linspace(min(y), max(y), options.ahatPts);
    p_ahat{i} = p/trapz(val, p);
    
    display(['Defect ' num2str(i)])
end

%% --- Find the threshold --- %
cum_cdf = cumtrapz(ahat{1}, p_ahat{1});
if isempty(options.Threshold) && isempty(options.PD)
    threshold = ahat{1}(find(cum_cdf > (1-options.FArate), 1));
    pf = [];
elseif ~isempty(options.PD)
    cum_cdf = cumtrapz(ahat{end}, p_ahat{end}); % given required PD for the biggest defect 
    threshold = ahat{end}(find(cum_cdf > (1-options.PD), 1));
    ind = find(ahat{1} > threshold, 1);
    pf = trapz(ahat{1}(ind:end), p_ahat{1}(ind:end));
else
    threshold = options.Threshold;
    ind = find(ahat{1} > threshold, 1);
    pf = trapz(ahat{1}(ind:end), p_ahat{1}(ind:end)); 
end

if isempty(threshold)
    threshold = ahat{1}(end);
end

pod = zeros(options.aPts-1, 1);
plow = zeros(options.aPts-1, 1);
phigh = zeros(options.aPts-1, 1);
for i=2:options.aPts
   cum_cdf = cumtrapz(ahat{i}, p_ahat{i});
   ind = find(ahat{i} > threshold, 1);
   pod(i-1) = trapz(ahat{i}(ind:end), p_ahat{i}(ind:end)); 
   plow(i-1) = pod(i-1) - 1/(4*options.samplePts*0.05^2);
   phigh(i-1) = pod(i-1) + 1/(4*options.samplePts*0.05^2);   
end

ai = ai(2:end);
