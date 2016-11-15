function pout = arrayVoting(pod, kn)

nsense = size(pod, 2);
ndefect = size(pod, 1);

% total = factorial(nsense)/(factorial(nsense-k)*factorial(k));
pout = zeros(ndefect,1);
% n = 1;

for i=1:ndefect
    for k=kn:nsense
        inds = 1:k;
        indsc = 1:nsense; 
        indsc(inds) = [];
        newInd = 2;
        compareInd = nsense-k+1:nsense;
        while  isempty(newInd) || newInd > 1
            pout(i) = pout(i) + prod(pod(i, inds)).*prod(1-pod(i,indsc));
%             display(['n=' num2str(n) ', ind= ' num2str(inds) ', indc= ' num2str(indsc)])
%              n=n+1;
            newInd = find((compareInd-inds) == 0,1);
            if isempty(newInd)
                inds(end) = inds(end)+1;                
            elseif newInd > 1
                tot = length(inds(newInd-1:end));
                inds(newInd-1:end) = inds(newInd-1)+1:inds(newInd-1)+tot;
            end
            indsc = 1:nsense;
            indsc(inds) = [];
        end
    end
end


