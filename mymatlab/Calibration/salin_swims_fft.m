function Sdata = salin_swims_fft(Cond, Temp, Pres, CStruct)
% usage: Sdata = salin_swims_fft(Cond, Temp, Pres, CStruct);
%  Cond [s/m], Temp [degC, in situ], Pres [dbar as of 2015],
%  CStruct.(lag,alfa,beta,freq, ??), computation parameters:
%     lag=scans to advance;
%     alfa,beta=thermal lag params;
%     freq = data rate [Def = 24 Hz];
%     ?? = future flags/values for filtering, etc.
%  Sdata.salin = salinity in PSU (e.g., 35.0);
%       .C_adj, .T_adj - filtered/adjusted C,T's for salin computation
% Dave Winkel, Jan-2005, 2015
% Changed to frequency-space shifting (May 2007, ashcherbina@apl.washington.edu)
clear Sdata
Sdata.salin = NaN * Pres; % error return values
%[size(Cond),size(Temp),size(Pres)]
if nargin>2 && length(Cond)==length(Temp) ...
        && length(Cond)==length(Pres) && length(Cond)>1
    x = 1;
else
    disp('salin_swims_fft: Need Cond,Temp,Pres vectors; returning NaNs')
    return
end
if nargin<4 || ~isstruct(CStruct)
    clear CStruct
    CStruct.lag = 0;
end

% defaults:
lag = 0; tau=0; alfa = 0; beta = 0; freq = 24;
vvs = {'lag','tau','alfa','beta','freq'};
% check for specified values
for i=1:length(vvs)
    if isfield(CStruct,vvs{i})
        x = CStruct.(vvs{i});
        if ~isempty(x)
            eval([vvs{i} ' = x;']);
        end
    end
end

% Shift, if able
if length(Cond)<5*freq,
    warning('Conductivity series too short?');
end

% convert legacy coefficients to the new ones (sign IS ok - DW)
lag   = lag/freq; % legacy lag is in scans -> convert to seconds
tau   =  tau; % there's no "tau" among legacy coeffs (ok as is)
alpha =  alfa;% alpha and beta coefficients are the same...
beta  =  beta;

% Salinity-spiking correction: match T and C sensor response
% 
if abs(lag)>1e-12 || abs(tau)>1e-12
    Cc = c_shift(Cond, lag, tau, 1/freq);
else
    Cc = Cond;
end

% Thermal mass correction: temperature used in salinity calculations should
% be filtered to simulate the effect of the conductivity cell
% heating(cooling) the water inside itself.
if abs(alpha)>1e-12 && abs(beta)>1e-12
    Tc = tmass_shift(Temp,alpha,beta,1/freq);
else
    Tc = Temp;
end

Sdata.salin = sw_salt(Cc*10/sw_c3515, Tc, Pres); % Temp s/b IPTS-68?
% jan-2012 - save corrected C,T's
Sdata.C_adj = Cc;
Sdata.T_adj = Tc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function C = c_shift(C,delta,tau,dt)
% perform filtering in freq. space
% [delta, tau]=[shift, time constant]
n = length(C);
C0 = mean(C); 
C = C-C0;
% pre-whitening (first difference)
C = diff(C');

[nfft,nt] = size(C);
C = fft(C);

nhf = floor(nfft/2);
if mod(nfft,2)
    % odd
    f=[0:nhf,-nhf:-1]'./(nfft*dt); 
    w = 2*pi*f;% radian
    F = exp(i*w*delta)./(1+i*w*tau); % filter function
else
    %even
    f=[0:nhf,-nhf+1:-1]'./(nfft*dt); 
    w = 2*pi*f;% radian
    F = exp(i*w*delta)./(1+i*w*tau); % filter function
    % fix Nyquist
    F(nhf+1) = real(F(nhf+1));
end

C = ifft(C.*repmat(F,1,nt),'symmetric');

C = [0; cumsum(C)];
C = C' - mean(C) + C0; % replace mean

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T = tmass_shift(T,alpha, beta ,dt)
% perform filtering in freq. space
% alpha: fractional "instanteneous" response 
% beta = 1/(time constant)

e = 1-alpha; 
tau = 1/beta;

% pre-whitening (first difference)
T0 = mean(T);
T = diff(T');

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

F = (1+e*iwt)./(1+iwt);% this converts measured temperature to that in the conductivity cell

T = T.*repmat(F,1,nt);
T = ifft(T,'symmetric');

T = [0; cumsum(T)];
T = T-mean(T);
T = T' + T0;
