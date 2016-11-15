function [peak, peaktime, tof, fcenter] = getSignalParams(t, y, tref, yref, varargin)

import Lamby.*

funName = 'getSignalParams';

% --- Default options values --- %
options = struct('filterCutoff', 0, ...
                 'filterType', 'lowpass', ...
                 'npts', length(t), ...
                 'sensitivity', 0.5);

% read the acceptable names
optionNames = fieldnames(options);

% count arguments
nArgs = length(varargin);
if round(nArgs/2) ~= nArgs/2
   error([funName ' needs propertyName/propertyValue pairs.'])
end

for pair = reshape(varargin,2,[]) %# pair is {propName;propValue}
   if any(strmatch(pair{1}, optionNames))
      % overwrite options. If you want you can test for the right class here
      % Also, if you find out that there is an option you keep getting wrong,
      % you can use "if strcmp(inpName,'problemOption'),testMore,end"-statements
      options.(pair{1}) = pair{2};
   else
      error('%s is not a recognized parameter name',pair{1})
   end
end


if options.filterCutoff ~= 0
    [t, y] = filterSignal(t, y, options.filterCutoff, options.npts,  options.filterType);
end

env = abs(hilbert(y, options.npts));
% env = abs(y);

[loc, pks] = peakfinder(env, max(env)*options.sensitivity);
% [pks, loc] = max(env);

if ~isempty(pks)
    peak = pks(1);
    peaktime = t(loc(1));
else
    peak = 0;
    peaktime = 0;
end

% reference peak 
env_ref = abs(hilbert(yref));
[loc_ref,  ~] = peakfinder(env_ref, max(env_ref)*options.sensitivity);

if peak < 1e-9
    tof = 0;
else
    tof = peaktime - tref(loc_ref(1));
end
    

% -- find the center frequency 
if length(pks) > 1
    sep = ceil((loc(2)+loc(1))/2);
    y(sep:end) = 0;
end

[f, Y] = findfft(t, y, options.npts, 'single');
df = (f(2) - f(1))/10;
fi = 0:df:f(end);
Y = interp1(f, abs(Y), fi, 'spline');

% figure
% plot(fi, angle(Yr+sqrt(-1)*Yi))

[~,indf] = max(Y);
fcenter = fi(indf);


end

