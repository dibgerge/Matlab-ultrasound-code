function [pod, threshold, plow, phigh] = findThreshold(ahat, p_ahat, pf, nsamples, confidence)

cum_cdf = cumtrapz(ahat{1}, p_ahat{1});
threshold = ahat{1}(find(cum_cdf > (1-pf), 1));
plow = [];
phigh = [];
if isempty(threshold)
    threshold = ahat{1}(end);
end

pod = zeros(length(ahat)-1, 1);
% plow = zeros(length(ahat}-1, 1);
% phigh = zeros(length(ahat}-1, 1);

for i=2:length(ahat)
   cum_cdf = cumtrapz(ahat{i}, p_ahat{i});
   ind = find(ahat{i} > threshold, 1);
   pod(i-1) = trapz(ahat{i}(ind:end), p_ahat{i}(ind:end));    
   
   if nargin > 3
        plow(i-1) = pod(i-1) - 1/(4*nsamples*(1-confidence)^2);
        phigh(i-1) = pod(i-1) + 1/(4*nsamples*(1-confidence)^2);   
   end
end

pod = pod(:);