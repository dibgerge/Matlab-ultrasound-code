function [x,t] = inversefouriertransform(Cx,varargin)  
%INVERSEFOURIERTRANSFORM   Performs the Real Inverse Fourier Transform.
%   [X,t] = INVERSEFOURIERTRANSFORM(RC,IC,dF) Gives the real time series 
%   X(t), from the components of the Fourier Complex Transform from zero to 
%   Nyquist frequencies  
%                CX(Ff) = RC(Ff) + i*IC(Ff) = COMPLEX(RC,IC)
%   with sampling interval dF, via Inverse Fast Fourier Transform.
%
%   [X,t] = INVERSEFOURIERTRANSFORM(CX,dF) Does the same thing. 
%
%   [X,t] = INVERSEFOURIERTRANSFORM(CX,dF,'exp') Gives the same thing but
%   via complex exponentials (not so fast as IFFT).
%
%   [X,t] = INVERSEFOURIERTRANSFORM(CX,dF,'trig') Gives the same thing via
%   cosine-sine series (not so fast as IFFT).
%
%   Note: because the transform of real series is conjugate even (the real 
%   part even and the imaginary part odd) the input is only for the 
%   positive Fourier frequencies:  
%                        Ff = (0:1/N:1/2)/dT  
%   i.e., from zero to the Nyquist frequency. Where N = 2*length(IC)-1 - 
%   (IC(end)==0) is the length of the time series, and dT = 1/(N*dF) his 
%   sampling interval.
%
%   Note: if the sampling interval is 0.5 cycles/minute, for example, and 
%   the X units, [X], are meters, then:  
%   a) PERCIVAL & WALDEN (normal usage)  
%       if  dF=1/120          =>  [RC,IC]=m*sec,   [t]=sec  
%   b) BLOOMFIELD (normalized time interval)  
%       if  dF=1              =>  [RC,IC]=m,       [t]=2min   
%   c) MATLAB (no dT input) (normalized frequency interval)  
%       if  dF=1/N or empty   =>  [RC,IC]=m*2min,  [t]=1 (adimentional)  
%     
%   Note: This program was created for teaching purposes. That's the reason 
%   of the 'exp' and 'trig' methods.  
%  
%   Example:  
%      N = 100; dT = 120; t = (0:N-1)*dT; f1 = 0.0007; f2 = 0.0013;  
%      x = 20*sin(2*pi*f1*t) + 30*cos(2*pi*f2*t) + rand(size(t));  
%      x = x-mean(x);  
%      [rc,ic,Ff] = fouriertransform(x,dT);
%      [xi,ti] = inversefouriertransform(rc,ic,1/(N*dT));
%      subplot(211), plot(t,x,'.-',ti,xi,'o'), 
%      xlabel('time, sec'), ylabel('x, m')  
%      legend('Original','Inverse')
%      title('Complex inverse Fourier transform example')
%      subplot(212), plot(Ff,rc,Ff,ic,[f1 f2],[0 0],'ro')  
%      xlabel('frequency, cycle/sec'), ylabel('Cx, m sec')  
%      legend('Real part','Imaginary Part','Natural Frequencies',4)
%  
%   See also FOURIERTRANSFORM, FFT, IFFT  
  
%   Written by  
%   Lic. on Physics Carlos Adrián Vargas Aguilera  
%   Physical Oceanography MS candidate  
%   UNIVERSIDAD DE GUADALAJARA   
%   Mexico, 2004  
%  
%   nubeobscura@hotmail.com  
  
% Time series information:  
[Cx,dF,N,N2,Method] = check_arguments(Cx,varargin,nargin);  
  
% Real inverse Fourier transform:
switch lower(Method)
 case 'fft'
  x = inversefouriertransform_fft(Cx,dF,N,N2);
 case 'exp'
  x = inversefouriertransform_exponential(Cx,dF,N,N2);
 case 'trig'
  x = inversefouriertransform_trigonometric(Cx,dF,N,N2); 
 otherwise
  error('Method unknown. Must be one of ''fft'', ''exp'' or ''trig''.')
end
  
% Time argument:  
t = (0:N-1)/(N*dF);     t = reshape(t,size(x));  

function x = inversefouriertransform_fft(Cx,dF,N,N2)  
% Real inverse Fourier transform via FFT.  
Cx(N2+1:N) = conj(Cx(N-N2+1:-1:2)); 
x = ifft(Cx)*N*dF;  x = real(x);    
  
function x = inversefouriertransform_exponential(Cm,dF,N,N2)  
% Real inverse Fourier transform via complex series.  
n  = 0:N-1;       m  =  n;  
tn = n;           % tn*dT:  time    
fm = m/N;         % fm/dT:  Fourier frequency  
wm = 2*pi*fm;     % angular Fourier frequency  
Cm(N2+1:N) = conj(Cm(N-N2+1:-1:2)); 
x = zeros(1,N);
for n = 1:length(tn)  
 x(n) = sum( Cm(:) .* exp(i*wm(:)*tn(n)) ) * dF ;  
end
x = real(x); [mc,nc] = size(Cm); x = shiftdim(x,mc>nc);

function x = inversefouriertransform_trigonometric(Cx,dF,N,N2)  
% Real inverse Fourier transform via cosine-sine series.
n  = 0:N-1;       m  =  0:N2-1;  
tn = n;           % tn*dT:  time    
fm = m/N;         % fm/dT:  Fourier frequency  
wm = 2*pi*fm;     % angular Fourier frequency  
% Translation from complex transform:
am =  2*real(Cx)*dF;
bm = -2*imag(Cx)*dF;
am(1) = am(1)/2;
if ~rem(N,2)
 am(end) = am(end)/2;
end
x = zeros(1,N);
for n = 1:length(tn)
 x(n) = sum( am(:).*cos(wm(:)*tn(n)) + bm(:).*sin(wm(:)*tn(n)) );
end
[mc,nc] = size(Cx); x = shiftdim(x,mc>nc);

function [Cx,dF,N,N2,Method] = check_arguments(Cx,Ventries,Nentries)  
% Is vector?  
N2 = length(Cx);  
if N2~=numel(Cx)  
 error('Entry must be a vector.')  
end  
% Sampling interval?  
ic = [];          % Default (complex entry)
dF = [];          % Default: 1/N (MATLAB way)
Method = 'fft';   % Default (MATLAB way)
dF2 = [];
if Nentries > 1
 if ischar(Ventries{1})
  Method = Ventries{1};
 elseif length(Ventries{1}) == N2
  ic = Ventries{1};
 else
  dF = Ventries{1};
 end
end
if Nentries > 2
 if ischar(Ventries{2})
  Method = Ventries{2};
 elseif length(Ventries{2}) == N2
  ic = Ventries{2};
 else
  dF = Ventries{2};
 end
end
if Nentries > 3
 if ischar(Ventries{3})
  Method = Ventries{3};
 elseif length(Ventries{3}) == N2
  ic = Ventries{3};
 else
  dF = Ventries{3};
 end
end
if ~isempty(ic)
 Cx = complex(Cx,ic);
end

% Series length (if the last imaginary element is zero, then is even):
N = 2*N2 - 1 - (imag(Cx(end))==0);

% Default sampling frequency interval:
if isempty(dF)
 dF = 1/N;
end


% Carlos Adrián. nubeobscura@hotmail.com