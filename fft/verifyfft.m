clear all 
close all

gamma = 2*pi*200;
a = 4e-3;
dx = 0.01e-3;
fs = 1/dx;
N = 40000;
x = (-(N-1)/2:(N-1)/2)*dx;

[C ind1] = min(abs(x+a));
[C ind2] = min(abs(x-a));

y = zeros(N, 1);
y(ind1:ind2) = sinh(gamma*x(ind1:ind2));
figure; plot(x, y);

dxi = 1/(x(end)-x(1)); 
xi = (-(N-1)/2:(N-1)/2)*dxi;
Ynum = dx*fftshift(fft(fftshift(y)));

Y = (sqrt(-1)*2./(xi.^2+gamma^2)).*(xi*sin(gamma*a).*cos(xi*a) - gamma*cosh(gamma*a).*sin(xi*a));
figure; plotyy(xi, imag(Y), xi, imag(Ynum))
%% 
T = 0:0.1e-6:50e-6;
alpha = 1.3062e4;
wbar = 2.0134e5;

% 'Original' Data
ipt = exp(-alpha*T).*sin(wbar*T);

% Take fft of ipt
Fs = 1/(T(2)-T(1));
L = length(ipt);
k = -(L-1)/2:(L-1)/2;

freq = k*Fs/L;

Ipw = fftshift(fft(ipt)/L);

% Fourier transform of 'original' data
w = 2*pi*freq;
Ip = wbar./((alpha+1i.*w).^2+wbar.^2);

% Compare both fourier transforms
figure
[AX, H1, H2] = plotyy(freq, Ipw, freq, Ip);
set(AX(2),'ylim',[-1.35e-5 1.35e-5])
set(AX(1),'xlim',[-1e5 1e5])
set(AX(2),'xlim',[-1e5 1e5])
set(get(AX(1),'Ylabel'),'String','Ipw')
set(get(AX(2),'Ylabel'),'String','Ip')

%% Verify a simple sine wave
fi = 10;
Fs = 5e2;
dt = 1/Fs;
cycles = 100;
T = cycles/fi;
N = T*Fs;
t = (0:N-1)*dt;
y = sin(2*pi*fi*t);

%shift back time axis, and append a zero at the end of signal
NFFT = length(t);
% NFFT = 16*2^nextpow2(L); % Next power of 2 from length of y
df = Fs/NFFT;
f = (0:(NFFT-1))*df;
f(f >= Fs/2) = f(f >= Fs/2) - Fs;
Y = fft(y, NFFT);
Y = fftshift(Y)/NFFT;
f = fftshift(f);

% figure; plot(t, y);
% hold on;
figure; plot(f, real(Y), f, imag(Y), f, abs(Y)) ;
% figure; plot(fftshift(f), fftshift(abs(Y)), '-*');
xlim([fi/2 fi*2]);

% plotyy(t*1e6, y, t2*1e6, y2);