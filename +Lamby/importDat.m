function [t, A0, S0] = importDat(filename, interp_factor, tmax)

data = load(filename);
t = data(:,1);
A0 = data(:,2:2:end) + data(:,3:2:end);
S0 = data(:,2:2:end) - data(:,3:2:end);

for i=1:length(A0(1,:))
    [~, A0(:,i)] = Lamby.filterSignal(t, A0(:,i), 300e3);
    [~, S0(:,i)] = Lamby.filterSignal(t, S0(:,i), 300e3);
end

if nargin < 3
    tmax = t(end);
end

if nargin > 1
    if interp_factor <= 0 || ~isscalar(interp_factor)
        error('interp_factor should be a positive scalar.');        
    end
    dt = t(2)-t(1);
    dti = dt/interp_factor; 
    N = floor(tmax/dti);
    ti = (0:N-1)*dti;
    
    A0 = interp1(t, A0, ti, 'spline');
    S0 = interp1(t, S0, ti, 'spline');
    t = ti;
end

end

