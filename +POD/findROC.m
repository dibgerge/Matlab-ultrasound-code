function [pf, pd] = findROC(ahat, p_ahat)

pf = zeros(length(ahat{1}),1);
pd = zeros(length(ahat{1}),length(ahat)-1);

for i=1:length(ahat{1})-1
    threshold = ahat{1}(i);
    pf(i) = trapz(ahat{1}(i:end).', p_ahat{1}(i:end).');    
    for j=2:length(ahat)
        ind = find(ahat{j} > threshold,1);
        if isempty(ind)
            pd(i,j-1) = 0;
        elseif ind >= length(ahat{j})
            pd(i,j-1) = 0;
        else
            pd(i,j-1) = trapz(ahat{j}(ind:end).', p_ahat{j}(ind:end).');
        end
    end
end


