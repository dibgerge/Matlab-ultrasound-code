%% 
% Inputs: y - The input signal
%         Threshold - the threshold to compute bandwidth or time window (dB)
%         Default value -30 dB
function [startInd, endInd] = findSignalLimit(y, threshold)

if ~isnumeric(threshold) || ~isscalar(threshold) || threshold >= 0
    error('Threshold should be a negative scalar (dB).');
end

if isvector(y) 
    y = [y(:); zeros(length(y),1)];
else
    y = [y; zeros(size(y))];
end

mult = 10^(threshold/10);
maxAmp = max(abs(y));

inds = find(y > maxAmp*mult);

if isempty(inds)
    startInd = [];
    endInd = [];
elseif length(inds) == 1
    startInd = inds(1);
    endInd = inds(1);
else    
    startInd = inds(1);
    endInd = inds(end);
end

end

