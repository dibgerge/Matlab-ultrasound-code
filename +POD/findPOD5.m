function [pod, plow, phigh, ai, ahat, p_ahat] = findPOD5(F, varargin)

% --- Default options values --- %
options = struct('PDFParams', [0 1], ...
                 'maxDefect', 10e-3, ...
                 'xiRange', [0 2], ...
                 'NoiseLevel', 0,...
                 'FArate', 0.05,...
                 'aPts', 100, ...
                 'ahatPts', 100,...
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
p_xi = makedist('Normal', 'mu', options.PDFParams(1), 'sigma', options.PDFParams(2));
p_xi = p_xi.truncate(options.xiRange(1), options.xiRange(end));

for i=1:options.aPts    
    %% --- (3) Get samples and find the 'ahat' PDF  ---%
    xivals = random(p_xi, options.samplePts, 1);
    
    [x, y] = ndgrid(xivals, ai(i));
    ahat_vals = F(x, y);
    
    [p, val] = hist(ahat_vals, options.ahatPts);    
    val = linspace(min(ahat_vals), max(ahat_vals), options.ahatPts);
    
    tot = trapz(val, p);
    p = p/tot;
    
    %% --- (4) Addd noise to the problem --- %
    if options.NoiseLevel > 0        
        nsigma = 5;
        noiseLimit = nsigma*options.NoiseLevel;
        dn = mean(diff(val));
        noiseRange = 0:dn:noiseLimit;
        
        noisePDF = makedist('Normal', 'mu', 0, 'sigma', options.NoiseLevel);
        p_xi = p_xi.truncate(0, noiseLimit);

        noise = pdf(noisePDF, noiseRange);
        convLen = length(val) + length(noiseRange) - 1;

        p_ahat{i} = conv(p(:), noise(:));

        ahat{i} = linspace(min(noiseRange)+min(val), ...
            max(noiseRange)+max(val), convLen); 

        p_ahat{i} = p_ahat{i}/trapz(ahat{i}, p_ahat{i});
    else
        p_ahat{i} = p;
        ahat{i} = val;
    end
    display(['Defect ' num2str(i)])
end

%% --- Find the threshold --- %
cum_cdf = cumtrapz(ahat{1}, p_ahat{1});
threshold = ahat{1}(find(cum_cdf > (1-options.FArate), 1));

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
