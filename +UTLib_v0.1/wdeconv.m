function G = wdeconv(Num,Den,NSR)
denom = (abs(Den)).^2; %norm(Den);
denom = denom+NSR;
denom = max(denom, sqrt(eps));
G = Num.*(conj(Den)./denom); %Correct 1 echo
return;
