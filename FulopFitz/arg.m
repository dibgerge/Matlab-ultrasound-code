% arg(X) where X is a matrix of complex numbers, computes the argument or phase angle of x
% but returns a value between 0 and 2pi, unlike the builtin Matlab angle function
% which returns the principle argument value as per standard between -pi and pi.

% the routine adds 2pi to negative angles returned by angle to accomplish this.
% the routine is needed for spectral analysis of matrices which have angular freq.
% encoded in compex phases; naturally the negative frequencies do not have the right
% physical interpretation.


function out = arg(X)

countj = 1;
for vec = X
    counti = 1;
    for x = vec.'
        princarg = angle(x);
 
        if sign(princarg) == -1
            out(counti,countj) = princarg + 2.*pi;
        else out(counti,countj) = princarg;
        end
        counti = counti + 1;
    end
    countj = countj + 1;
end

