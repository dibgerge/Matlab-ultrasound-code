function [t, y, f, Y] = hannpulse(fc, ncycle, varargin)
%HANNPULSE generate a hanning modulated tone burst
%   Detailed explanation goes here

% --- Default options values --- %
options = struct('Fs', fc*500, ...
                 'T', 0, ...
                 'writeFile', false, ...
                 'Npts', 0);

% read the acceptable names
optionNames = fieldnames(options);
% count arguments
nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2
   error('hannpulse needs propertyName/propertyValue pairs.')
end

for pair = reshape(varargin,2,[])
   if any(strmatch(pair{1}, optionNames))
      options.(pair{1}) = pair{2};
   else
      error('%s is not a recognized parameter name',pair{1})
   end
end

Ti = ncycle/fc;
N = Ti*options.Fs;

t = ((0:(N-1))/options.Fs)';
y = 0.5*(1-cos(2*pi*fc*t/ncycle)).*cos(2*pi*fc*t);

if options.T == 0 && options.Npts > 0
    options.T = options.Npts/options.Fs;
end

if Ti < options.T
    Nm = ceil((options.T - Ti)*options.Fs);
    t = [t; (N:(N+Nm-1)).'/options.Fs];
    y = [y; zeros(Nm, 1)];
end

[f, Y] = SPlib.findfft(y, options.Fs);

if options.writeFile
    vals = [t y];
    save(['hann_' num2str(fc/1e3) 'k_' num2str(ncycle) 'cyc.txt'], 'vals', '-ascii', '-tabs');
end
