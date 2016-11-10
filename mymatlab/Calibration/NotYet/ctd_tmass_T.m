function T = ctd_tmass_T(T,ab,dt)
% this should be applied to T

% pad with zeros (if fewer parameters passed)
ab = [ab,0 0 0];
ab = ab(1:3);

a = ab(1);
e = 1-a;
tau = ab(2);
delta = ab(3);

if nargin<3,
    dt = 1;
end

% pre-whitening (first difference)
T0 = mean(T);
T = diff(T);

[nfft,nt] = size(T);
T = fft(T);
% nhf = size(T,1);
nhf = floor(nfft/2);
if mod(nfft,2)
    % odd
    f=[0:nhf,-nhf:-1]'./(nfft*dt); % 1/s
else
    %even
    f=[0:nhf,-nhf+1:-1]'./(nfft*dt); %
end
w = 2*pi*f;% radian


iwt = (i*tau)*w;
if delta~=0,
    F = exp((-i*delta)*w).*(1+e*iwt)./(1+iwt);% this converts measured temperature to that in the conductivity cell
else
    F = (1+e*iwt)./(1+iwt);% this converts measured temperature to that in the conductivity cell
end
T = T.*repmat(F,1,nt);
T = ifft(T,'symmetric');

T = [0;cumsum(T)];
T = T-mean(T);
T = T+T0;
