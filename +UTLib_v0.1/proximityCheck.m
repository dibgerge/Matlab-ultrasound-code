function inprox = proximityCheck(pts, dist)
%assume 'pts' is monotonic vector

inprox = {};
i = 1;
n = 1;
while i < length(pts)
    d = abs(pts(i:end)  - pts(i));
    ind = find(d < dist);
    ind = reshape(ind, 1, []);
    
    if length(ind) > 1
        inprox{n} = i-1+ind; 
        i = ind(end)+i;    
        n = n+1;
    else
        i = i + 1;
    end        
end
end

