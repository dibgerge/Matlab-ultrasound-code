function [ahat, p_ahat, snr, mu, sigma] = findDistributionSingle(F, varargin)

% --- Default options values --- %
options = struct('PDFParams', [0 1], ...
                 'defectRange', 0:10, ...
                 'sampleBins', 100, ...
                 'xiRange', [0.1 2], ...
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
 
p_ahat = cell(length(options.defectRange),1);
ahat = cell(length(options.defectRange),1);

%% --- Obtain the variability PDF ---%
if options.PDFParams(2) > 0
    p_xi = makedist('Normal', 'mu', options.PDFParams(1), 'sigma', options.PDFParams(2));
    p_xi = p_xi.truncate(options.xiRange(1), options.xiRange(end));
end

mu = [];
sigma = [];
snr = cell(length(options.defectRange),1);

for i=1:length(options.defectRange)
    xivals = random(p_xi, options.samplePts, 1);    
    [xi, a] = ndgrid(xivals, options.defectRange(i));
    z = F(xi, a);
    [p, val] = hist(z, options.sampleBins);             
    ahat{i} = linspace(min(val), max(val), options.sampleBins);
    p_ahat{i} = interp1(val, p, ahat{i}, 'spline', 'extrap');
    p_ahat{i} = p_ahat{i}/trapz(ahat{i}, p_ahat{i});
        
    if options.defectRange(i) == 0
        mu = trapz(ahat{i}, ahat{i}.*p_ahat{i});
        sigma = sqrt(trapz(ahat{i}, p_ahat{i}.*(ahat{i}-mu).^2));
    end
        
    if ~isempty(mu)
        snr{i} = 20*log10(ahat{i}/sigma);
    end
    
    display(['Defect ' num2str(i)])        
end
