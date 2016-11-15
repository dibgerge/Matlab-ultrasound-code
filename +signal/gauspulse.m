% F: Gaussian pulse center frequency
%bw: Bandwidth limit where after it, it has 40dB attenuation
%    this is relative bandwidth: typical value 0.3
%outputFile: Select if an output file is needed with the signal info
%showplot (optional): select if time and fft plots will be drawn
function [t, y, f, Y] = gauspulse(F, bw, varargin)

import SPlib.*

% --- Default options values --- %
options = struct('dt', 1e-9, ...
                 'T', 0, ...
                 'threshold', -50, ...
                 'writefile', false,...
                 'showplot', false);

%--- read the acceptable names
optionNames = fieldnames(options);

%--- count arguments
nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2
   error('gauspulse needs propertyName/propertyValue pairs.')
end

for pair = reshape(varargin,2,[]) %# pair is {propName;propValue}
   if any(strncmp(pair{1}, optionNames, length(pair{1})))
      options.(pair{1}) = pair{2};
   else
      error('%s is not a recognized parameter name',pair{1})
   end
end

tc = gauspuls('cutoff', F, bw, [], options.threshold); 
t = (-tc : options.dt : tc)'; 
y = gauspuls(t, F, bw);

%shift back time axis
t = t+tc;

if  options.T > t(end)
    t2 = (0:options.dt:options.T)';
    y = interp1(t,y, t2, 'linear', 0);
    t = t2;
end

[f, Y] = findfft(t, y);
NFFT = length(Y);

if options.writefile
    vals = [t' y'];
    save(['gauss_' num2str(F/1e3) 'k_bw_' num2str(bw*100) '.txt'], 'vals', '-ascii', '-tabs');
end

%% plotting
if options.showplot         
    f = fftshift(f);
    Y = fftshift(Y);
    figure;
    plot(t*1e6, y);
    title(['$f(t) = cos(2\pi f_0 t)e$-${\frac{\pi^2 ln(10) (' num2str(bw)...
        'f_0)^2t^2}{5}}$, $f_0=' num2str(F/1e3) 'kHz$'], ...
        'Interpreter', 'latex',...
        'Units', 'normalized',...
        'Position', [0.5 0.99], ...
        'FontSize', 36)
    xlabel('\mus');
    ylabel('f(t)');
    xlim([0 (tc*2+0.4*tc)*1e6])
    ylim([-1.75 1.75])
    legend('Tone burst', 'H(\omega)', 'F(\omega)')

    figure;
    plot(f*1e-3, abs(Y)/NFFT, '-b');
%     title('$\tilde{f}( \omega )$', 'Interpreter', 'latex', 'Units', 'normalized', ...
%         'Position', [0.5 0.99],...
%         'FontSize', 36);
     xlabel('kHz');
%     ylabel('$|\tilde{f}(\omega)|$', 'Interpreter', 'latex');
     xlim([0 F*2*1e-3]);
end
