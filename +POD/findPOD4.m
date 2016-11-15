function [ahat, p_ahat] = findPOD4(xi, a, feature, varargin)

% --- Default options values --- %
options = struct('PDFParams', [0 1], ...
                 'defectRange', [0 10e-3], ...
                 'NoiseLevel', 0,...
                 'FArate', 0.05,...
                 'aPts', 100, ...
                 'ahatPts', 100);

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
 
%% --- Initialization and interpolate original axis --- %
ai = linspace(a(1), a(end), options.aPts);
xii = linspace(xi(1), xi(end), options.ahatPts);
p_ahat = zeros(options.ahatPts, options.aPts);
ahat = zeros(options.ahatPts, options.aPts);


%% --- (1) Construct the interpolation function --- % 
[X, Y] = ndgrid(xi, a);
F = griddedInterpolant(X, Y, feature');

%% --- (2) Obtain the variability PDF ---%
p_xi = makedist('Normal', 'mu', options.PDFParams(1), 'sigma', options.PDFParams(2));
p_xi = p_xi.truncate(xi(1), xi(end));


%%
for i=1:1%options.aPts
    % ---  Find the range of amplitudes due to the variability ---%
    [x, y] = ndgrid(xii, ai(i));
    amps  = F(x, y);
    
    % --- 
%     st = options.PDFParams(1) - 2*options.PDFParams(2);
%     en = options.PDFParams(1) + 2*options.PDFParams(2);
%     [x,y] = ndgrid(linspace(st,en, ceil(options.ahatPts/2)), ai(i)); 
%     a1 = F(x,y);
%     ahat(:,i) = linspace(min(a1), max(a1), options.ahatPts);
  
     ahat(:,i) = linspace(min(amps), max(amps), options.ahatPts);

    % --- Find changing direction point --- %
    n = 1;
    currStart = 1;
    
    for j=2:length(amps)-1
        if (amps(j) < amps(j-1) && amps(j) < amps(j+1)) || ...
                (amps(j) > amps(j-1) && amps(j) > amps(j+1))            
        % --- Construct the domains --- %
            domains(n) = makeDomain(currStart, j, amps, xii);
            currStart = j+1;            
            
        % --- Make the inverse function interpolant --- %
            Finv{n} = griddedInterpolant(domains(n).ahat, domains(n).xi);
            
        %  --- Find the derivative of the inverse functions --- %
            dxi = diff(domains(n).xi)./diff(domains(n).ahat);
            dahat = domains(n).ahat(1:end-1) + 0.5*(domains(n).ahat(2) - domains(n).ahat(1));
            dFinv{n} = griddedInterpolant(dahat, dxi, 'spline');
            
            n = n + 1;            
        end
    end
    
    domains(n) = makeDomain(currStart, length(amps), amps, xii);
    Finv{n} = griddedInterpolant(domains(n).ahat, domains(n).xi);    
    dxi = diff(domains(n).xi)./diff(domains(n).ahat);
    dahat = domains(n).ahat(1:end-1) + 0.5*(domains(n).ahat(2) - domains(n).ahat(1));
%     dahat = domains(n).ahat(1:end-1);
    dFinv{n} = griddedInterpolant(dahat, dxi, 'spline');
    
    for j=1:length(ahat(:,i))
        deriv_solns = [];
        func_solns = [];
        for k=1:length(domains)
            if ahat(j,i) <= domains(k).max && ahat(j,i) >= domains(k).min
                deriv_solns = [deriv_solns,  dFinv{k}(ahat(j,i))];
                func_solns = [func_solns,  Finv{k}(ahat(j,i))];
            end
        end
        p_ahat(j,i) = sum(abs(deriv_solns).*pdf(p_xi, func_solns));
    end

%     % --- (6) Compute the pdf function --- %
%         p_ahat(j,i) = sum(abs(1./dF(r)).*pdf(p_xi, r));
%     end
    display(['Defect ' num2str(i)])
end