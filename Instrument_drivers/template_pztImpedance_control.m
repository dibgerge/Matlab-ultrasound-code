%% Start the impedance analyzer
imp = HP4194();
imp.start();
imp.func = 1;
imp.impedance = 2;

%% Set-up for low frequency 
imp.averaging(16);
df1 = (20e3-1e3)/400;
f1 = 1e3:df1:20e3;
imp.sweep(1, [1e3, 20e3, 401]);
imp.scaleR([0, 1e3])
imp.scaleI([-50e3, 0])

%% Acquire low frequency
R1 = imp.getR();
C1 = imp.getI();
Z1 = R1 + sqrt(-1)*C1;

plot(f1, imag(Z1))

%% Set-up for high frequency
imp.averaging(8)
fstart = 20e3;
fend = 1e6;
df2 = (fend-fstart)/400;
f2 = fstart:df2:fend;
imp.sweep(1, [fstart, fend, 401]);
imp.scaleR([0, 2e3])
imp.scaleI([-4e3, 1.5e3])

%% Acquire high frequency
R2 = imp.getR();
C2 = imp.getI();
Z2 = R2 + sqrt(-1)*C2;

%%
imp.close();