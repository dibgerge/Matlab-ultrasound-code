function domain = makeDomain(startInd, endInd, ahat, xi)

domain.startInd = startInd;
domain.endInd   = endInd;
domain.min = min(ahat(startInd:endInd));
domain.max = max(ahat(startInd:endInd));
domain.ahat = ahat(startInd:endInd);
domain.ahat = domain.ahat(:);
domain.xi = xi(startInd:endInd);
domain.xi = domain.xi(:);

if domain.ahat(1) > domain.ahat(end)
    domain.ahat = flipud(domain.ahat);
    domain.xi = flipud(domain.xi);
end


            


end

