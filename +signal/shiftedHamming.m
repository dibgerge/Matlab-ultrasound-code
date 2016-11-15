function wind = shiftedHamming(len, centloc)
%SHIFTWINDOW shift a Hamming windowing function by the specified amount
%   - len: length of the hamming window
%   - cenloc: desired location of the window center

nshift = ceil(len/2)-centloc;
wind = hamming(len+2*abs(nshift));
if nshift > 0 
    wind = wind((2*nshift+1):end);
else
    wind = wind(1:(end-2*abs(nshift)));
end


end

