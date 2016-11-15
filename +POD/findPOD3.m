function [ahat, p_ahat] = findPOD3(xi, a, feature, varargin)

% --- Default options values --- %
options = struct('PDFParams', [0 1], ...
                 'defectRange', [0 10e-3], ...
                 'NoiseLevel', 0,...
                 'FArate', 0.05,...
                 'xiPts', 100,...
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
 
%% --- Initialization --- %
ai = linspace(a(1), a(end), options.aPts);
xii = linspace(xi(1), xi(end), options.xiPts);
p_ahat = zeros(options.xiPts, options.aPts);
ahat = zeros(options.ahatPts, options.aPts);


%% --- (1) Construct the interpolation function --- % 
[X, Y] = ndgrid(xi, a);
F = griddedInterpolant(X, Y, feature', 'spline');

%% --- (2) Obtain the variability PDF ---%
p_xi = makedist('Normal', 'mu', options.PDFParams(1), 'sigma', options.PDFParams(2));
p_xi = p_xi.truncate(xi(1), xi(end));


%%
for i=1:1
    % ---  (3) Find the range of amplitudes due to the variability ---%
    [x, y] = ndgrid(xii, ai(i));
    amps  = F(x, y);    
%     ahat(:,i) = linspace(min(amps), max(amps), options.ahatPts);
    [ahat(:,i), ind] = sort(amps);
    % --- (4) Find the derivate of the curve at a given defect size ---% 
    dvals = diff(amps)./diff(x);
    dx = x(1:end-1) + (x(2)-x(1))/2;
    dF = griddedInterpolant(x(1:end-1), dvals, 'spline');
 
    p_ahat(:,i) = abs(1./dF(xii(ind))).*pdf(p_xi, xii(ind));
%     ops = optimset('Display', 'off');
%     for j=1:options.ahatPts        
%     % --- (5) Find the solutions a_hat = f(xi,a) --- %
%         r = [];
%         fun = @(x)(F(x, ai(i))-ahat(j,i));    
%         for k=1:length(xi)
%             r = [r, fzero(fun, xi(k), ops)];
%         end
%         r = unique(r);
%         DUPE = diff([r NaN]) < 1e-12;
%         r(DUPE) = [];
%         r(isnan(r)) = [];
%         
%         if find(diff(r)) < 1e-6
%             error('one root');
%         end
%             
%     % --- (6) Compute the pdf function --- %
%         p_ahat(j,i) = sum(abs(1./dF(r)).*pdf(p_xi, r));
%     end
    display(['Defect ' num2str(i)])
end