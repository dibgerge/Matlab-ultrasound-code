function y = removePulse(t, pulse, response)
y =[];
%%________________________________________________________________________%

import Lamby.*

[fp, Yp] = findfft(t, pulse);
[fr, Yr] = findfft(t, response);

ind1 = find(fp > 100e6, 1);
ind2 = find(fp < -100e6,1,'last');

Y = [Yr(1:ind1)./Yp(1:ind1); zeros(ind2-ind1-1,1); ...
    Yr(ind2:end)./Yp(ind2:end)];

%plot(fp, abs(Y))

y = ifft(Y);


%xlim([0 200e6])

