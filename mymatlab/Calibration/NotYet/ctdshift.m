function C = ctdshift(C,ab,dt)
% CTDSHIFT - apply sensor correction to a CTD time series
%   Cnew = ctdshift(C,ab,dt) applies the correction given by ab = [delta, tau]=[shift, time constant]
%   The transformation is equivalent to 
%       1. Advancing C by (delta) seconds 
%       2. Applying exponential filter with a time constant (tau) seconds
%   To perform salinity de-spiking on a typical CTD, ctdshift should be applied to conductivity with both coefficients being positive
%   (typical values for 24Hz data with no deck box shifting are [0.08 0.031] seconds or [1.9199  0.9152] scans)

ab = [ab(:);0;0];
ab = ab(1:2);
if nargin<3
    % assume unit timestep
    dt = 1;
end

% FFT-iFFT method (convolution is done through fft, anyway!)
n = length(C);
C0 = mean(C); 
C = C-C0;

% pre-whitening (first difference)
C = diff(C);

[nfft,nt] = size(C);
C = fft(C);

nhf = floor(nfft/2);
if mod(nfft,2)
    % odd
    f=[0:nhf,-nhf:-1]'./(nfft*dt); 
    w = 2*pi*f;% radian
    F = exp(i*w*ab(1))./(1+i*w*ab(2));
else
    %even
    f=[0:nhf,-nhf+1:-1]'./(nfft*dt); 
    w = 2*pi*f;% radian
    F = exp(i*w*ab(1))./(1+i*w*ab(2));
    % fix Nyquist
    F(nhf+1) = real(F(nhf+1));
end


C = ifft(C.*repmat(F,1,nt),'symmetric');

C = [0;cumsum(C)];
C = C-mean(C)+C0; % replace mean