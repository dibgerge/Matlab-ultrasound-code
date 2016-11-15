function [F, X, Y] = makeFunction(xi, a, feature)

%% --- Construct the interpolation function --- % 
[X, Y] = ndgrid(xi, a);
F = griddedInterpolant(X, Y, feature, 'spline');


end

