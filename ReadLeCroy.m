function u = ReadLeCroy(fname)
%READLECROY Summary of this function goes here
%   Detailed explanation goes here
w = ReadLeCroyBinaryWaveform(fname);
u = UTlib.utsignal(w.x, w.y);

end

