function [startInd, tstart] = findTriggerTime(t, y, fc, BW)

import Lamby.*
[ta, ya, ~, ~] = gauspulse(fc, BW, [], 0, 0);
[~, Ia] = max(ya);
ta_peak = ta(Ia);

[~, I] = max(y);
t_peak = t(I);

tstart = t_peak - ta_peak;

[~, startInd] = min(abs(t - tstart));



