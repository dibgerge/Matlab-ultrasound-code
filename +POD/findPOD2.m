function [ahat, p_ahat] = findPOD2(f, varargin)

% --- Default options values --- %
options = struct('PDFParams', [0 1], ...
                 'VariabilityRange', [0 1], ...
                 'defectRange', [0 10e-3], ...
                 'NoiseLevel', 0,...
                 'FArate', 0.05);

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

% --- Initialization ---%
defects = linspace(options.defectRange(1), options.defectRange(2));
variability = linspace(options.VariabilityRange(1), options.VariabilityRange(2));
p_ahat = zeros(100);
ahat = zeros(100);

% --- (1) Obtain the variability PDF ---%
pxi = makedist('Normal', 'mu', options.PDFParams(1), 'sigma', options.PDFParams(2));
pxi = pxi.truncate(options.VariabilityRange(1), options.VariabilityRange(2));

% --- (2) Find the derivative of the polynomial, assuming 4th order polynomial...
c = coeffvalues(f);
domains = [];
for k=1:length(defects)
    a = defects(k);
    
    pv = [c(11), ...
            c(7) + c(12)*a, ...
            c(4) + c(8)*a + c(13)*a^2, ...
            c(2) + c(5)*a + c(9)*a^2 + c(14)*a^3, ...
            c(1) + c(3)*a + c(6)*a^2 + c(10)*a^3 + c(15)*a^4];
    
    dpoly = [4*c(11), ...
            3*(c(7) + c(12)*a), ...        
            2*(c(4) + c(8)*a + c(13)*a^2), ...
            c(2) + c(5)*a + c(9)*a^2 + c(14)*a^3];

    d2poly = [12*c(11), ...
            6*(c(7) + c(12)*a), ...        
            2*(c(4) + c(8)*a + c(13)*a^2)];
        
    % Find the real roots and monotonous domains
%     rts = roots(dpoly);
%     for i=1:length(rts)
%         if isreal(rts(i)) && polyval(d2poly, rts(i)) ~= 0
%             domains = [domains, rts(i)];
%         end
%     end
    
    % find the limits of a_hat under the limits of the variability
%     festimate = feval(f, variability, repmat(a, 1, length(variability)));
    festimate  = polyval(pv, variability);
    fmin = min(festimate)*1.005;
    fmax = max(festimate);
    ahat(:,k) = linspace(fmin, fmax);
    
    %
    for i=1:length(ahat(:,k))
        pv_temp = pv;
        
        % find the real solutions for the equation ahat = f(xi)
        pv_temp(end) = pv_temp(end) - ahat(i,k);
        rts = roots(pv_temp);
        rts = rts(imag(rts) == 0);

        
        % find the ahat pdf at that point
        p_ahat(i,k) = sum(abs(1./polyval(dpoly,rts)).*pdf(pxi, rts));

    end

end

