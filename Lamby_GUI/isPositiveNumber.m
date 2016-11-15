function isvalid = isPositiveNumber(num)

isvalid = 0;

if isscalar(num) && isnumeric(num) && num >= 0
    isvalid = 1;
end

end

