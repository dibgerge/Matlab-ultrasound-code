function dispersionRemoval(t, sig, f0)
import Lamby.*

Fs = 1/(t(2)-t(1));
NFFT = length(t);
% NFFT = 8*2^nextpow2(L); % Next power of 2 from length of y
df = Fs/NFFT;
f = (0:(NFFT-1))*df;
f(f >= Fs/2) = f(f >= Fs/2) - Fs;
Y = fft(sig,NFFT);

%get only the single sided fft for now
indright = find(f>= 0);           
fright = f(indright);
Yright = Y(indright);
% compute only the response due to components with value larger
% than 0.1%
inds = find(abs(Yright) > 1e-3*max(abs(Yright)));

% [cs fInds] = plate.symmetricPhaseVelocity(fright(inds), 0.05);
[ca fInda] = plate.antisymmetricPhaseVelocity(fright(inds), 0.05);
% [cgs fIndgs] = plate.symmetricGroupVelocity(fright(inds), 0.01);
[cga fIndga] = plate.antisymmetricGroupVelocity(fright(inds), 0.01);

ka = 2*pi*fright(inds)./ca{1};
[C ind0] = min(abs(fright(inds) - f0));
kaLin =2*pi*fdraw(ind0)/ca{1}(ind0) + (2*pi*fdraw - 2*pi*f0)./cga{1}(ind0);

figure; plot(f(inds), ka, f(inds), kaLin)
figure; plot(t, 
end

