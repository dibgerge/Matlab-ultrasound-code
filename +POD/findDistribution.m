function [ahat, p_ahat] = findDistribution(F, varargin)

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

for i=1:length(options.defectRange)
    xivals = random(p_xi, options.samplePts, 1);    
    [xi, a] = ndgrid(xivals, options.defectRange(i));    
    
    p = zeros(options.sampleBins, length(F))';
    val = zeros(options.sampleBins, length(F))';
    minval = 100;
    maxval = 0;
    for j=1:length(F)
        z = F{j}(xi, a);
        [p(j,:), val(j,:)] = hist(z, options.sampleBins); 
        
        if min(val(j,:)) < minval
            minval = min(val(j,:));
        end
        if max(val(j,:)) > maxval
            maxval = max(val(j,:));
        end        
%         Fp = griddedInterpolant(val, p, 'spline', 'none');
%         p = Fp(options.featureRange);
%         p(isnan(p)) = 0;
    end
    
    ahat{i} = linspace(minval, maxval, options.sampleBins);        
    p_ahat{i} = ones(options.sampleBins,1);
    for j=1:length(F)
        p(j,:) = interp1(val(j,:), p(j,:), ahat{i}, 'linear', 0);
%         if i==1
%         plot(ahat{i}, p(j,:)/trapz(ahat{i}, p(j,:)));
%         hold on;
%         end
        p_ahat{i} = p_ahat{i}.*p(j,:)'; 
        p_ahat{i} = p_ahat{i}/trapz(ahat{i}, p_ahat{i});        
    end
    
    
    display(['Defect ' num2str(i)])
end
