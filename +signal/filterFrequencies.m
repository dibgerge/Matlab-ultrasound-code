function Y = filterFrequencies(f, Y, frange, option)

if nargin == 3
    option = 'lowpass';
end

switch(option)
    case 'lowpass'
        if ~isscalar(frange) || frange <= 0 
            error('frange should be a scalar and larger than zero.');
        end
        indrange = [find(f > frange); find(f < -frange)];
    case 'highpass' 
        if ~isscalar(frange) || frange <= 0 
            error('frange should be a scalar and larger than zero.');
        end
        indrange = find(f < frange & f > -frange);
    case 'bandpass'
        if length(frange) ~= 2 ||  frange(2) <= frange(1)
            error('frange should have two positive components, that are larger than zero.');
        end
        indrange = [find(f >= frange(2)); find(f <= -frange(2))];
        indrange = [indrange; find(f <= frange(1) & f >= -frange(1))];
    case 'notch'
         if ~isscalar(frange) || frange <= 0 
             error('frange should be a scalar and larger than zero.');
         end
    otherwise
        error('Unknown option.');
end
        
            
Y(indrange) = 0;


end

