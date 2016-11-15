function [feature_out, featurePDF, variabilityPDF, pod] = ...
    findPOD( variability, feature, varargin)

% supportedFits = {'poly', 'logpoly'};

% --- Default options values --- %
options = struct('Fit', 'logpoly', ...
                 'PolyOrder', 3, ...
                 'PDF', 'Normal', ...
                 'PDFParams', [0 max(variability)/2], ...
                 'NoiseLevel', 0,...
                 'Threshold', 0.05);

%# read the acceptable names
optionNames = fieldnames(options);

%# count arguments
nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2
   error('findPOD needs propertyName/propertyValue pairs.')
end

for pair = reshape(varargin,2,[]) %# pair is {propName;propValue}
   if any(strmatch(pair{1}, optionNames))
      %# overwrite options. If you want you can test for the right class here
      %# Also, if you find out that there is an option you keep getting wrong,
      %# you can use "if strcmp(inpName,'problemOption'),testMore,end"-statements
      options.(pair{1}) = pair{2};
   else
      error('%s is not a recognized parameter name',pair{1})
   end
end

% ------------------------------------------------------------------------%
% ------------------------------------------------------------------------%

% --- Compute the variability PDF function --- %
if length(options.PDFParams) == 1
    variabilityPDF = pdf(options.PDF, variability, options.PDFParams(1));
elseif length(options.PDFParams) == 2
    variabilityPDF = pdf(options.PDF, variability, options.PDFParams(1), ...
        options.PDFParams(2));
elseif length(options.PDFParams) == 3
    variabilityPDF = pdf(options.PDF, variability, options.PDFParams(1), ...
        options.PDFParams(2), options.PDFParams(3));
else
    error('PDFParameters has too many elements.');
end


% -- set variability parameter to its log if logpoly fit is required
if strcmpi(options.Fit, 'poly')
    [poly_coefs, S] = polyfit(variability, feature, options.PolyOrder);

    % evaluate the actual fitting function
    fittedFeature = polyval(poly_coefs, variability);

    % get the derivate of the inverse of the fitted polynomial function
    powers = fliplr(1:options.PolyOrder);
    inverseDeriv  = 0;
    for i=1:options.PolyOrder
        inverseDeriv = inverseDeriv + powers(i)*poly_coefs(i)*variability.^(powers(i)-1);
    end
    inverseDeriv = abs(1./inverseDeriv);
% - log polynomial fitting - %    
elseif strcmpi(options.Fit, 'logpoly')
    
    [poly_coefs, S] = polyfit(log(variability), feature, options.PolyOrder);

    % evaluate the actual fitting function
    fittedFeature = polyval(poly_coefs, log(variability));

    % get the derivate of the inverse of the fitted polynomial function
    denominator = 0;
    for i=1:options.PolyOrder
        denominator = denominator + (options.PolyOrder-i+1)*poly_coefs(i)...
            *log(variability).^(options.PolyOrder-i);
    end
    inverseDeriv = abs(variability./denominator);
end                
    
% multiply the derivative with the pdf function to compute the feature PDF
featurePDF = inverseDeriv.*variabilityPDF;
% feature_out = fittedFeature;
%interpolate to make the amplitude vector uniform 
feature_out = linspace(min(fittedFeature), max(fittedFeature), length(featurePDF));
featurePDF = interp1(fittedFeature, featurePDF, feature_out);

% --- Compute the noise PDF function --- %
if options.NoiseLevel > 0
    nsigma = 5;
    noiseLimit = nsigma*options.NoiseLevel;
    nnoise = 100;    
    dn = 2*noiseLimit/nnoise;
        
    if mean(diff(feature_out)) > dn
        feature_out_temp = min(fittedFeature):dn:max(fittedFeature);
        featurePDF = interp1(feature_out, featurePDF, feature_out_temp);
        feature_out = feature_out_temp;
    else
        dn = mean(diff(feature_out));
    end

    noiseRange = -nsigma*options.NoiseLevel:dn:nsigma*options.NoiseLevel;
    noisePDF = normpdf(noiseRange, 0, options.NoiseLevel);
    
    convLen = length(feature_out) + length(noiseRange) - 1;
    
    featurePDF = conv(featurePDF(:), noisePDF(:));
 
    feature_out = linspace(min(noiseRange)+min(feature_out), ...
        max(noiseRange)+max(feature_out), convLen); 
    
    normalization_const = trapz(feature_out, featurePDF);

    featurePDF = featurePDF/normalization_const;
end

startind = find(feature_out > options.Threshold, 1);
pod = trapz(feature_out(startind:end), featurePDF(startind:end));

end

