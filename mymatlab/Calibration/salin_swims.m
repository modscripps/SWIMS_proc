function Sdata = salin_swims(Cond, Temp, Pres, CStruct)
% usage: Sdata = salin_swims(Cond, Temp, Pres, CStruct)
%  Cond [s/m], Temp [degC, in situ], Pres [dbar - 6/2015],
%  CStruct.(lag,alfa,beta,freq, ??), computation parameters:
%     lag=scans to advance;
%     alfa,beta=thermal lag params;
%     freq = data rate [Def = 24 Hz];
%     ?? = future flags/values for filtering, etc.
%  Sdata.salin = salinity in PSU (e.g., 35.0);
%     FUTURE: return filtered/adjusted cond,temp records: .C_adj,.T_adj
% Dave Winkel, Jan-2005, 2015

clear Sdata
Sdata.salin = NaN * Pres; % error return values
%[size(Cond),size(Temp),size(Pres)]
if nargin>2 && length(Cond)==length(Temp) ...
        && length(Cond)==length(Pres) && length(Cond)>1
    x = 1;
else
    disp('salin_swims: Need Cond,Temp,Pres vectors; returning NaNs')
    return
end
if nargin<4 || ~isstruct(CStruct)
    clear CStruct
    CStruct.lag = 0;
end

% defaults:
lag = 0; alfa = 0; beta = 0; freq = 24;
vvs = {'lag','alfa','beta','freq'};
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
Craw = Cond;
if length(Cond) > abs(lag)+2
    Cond = shift(Craw, lag);
end

%% Thermal Mass algorithm is from SeaBird SeaSoft-Win32 manual (pdf)
if alfa>1e-10 && beta>=0
    % compute/initialize temp diffs, cond corrections
    dTp = Temp;
    dTp(2:end) = diff(Temp);
    dTp(1) = dTp(2);
    dcdt = 0.1 * (1 + 0.006*(Temp-20));
    ctm = 0*dTp;
    % a,b
    aa = 2 * alfa / (2 + beta/freq);
    bb = 1 - (2*aa/alfa);
    % compute corrections
    for i=2:length(Cond)
        ctm(i) = -1.0*bb*ctm(i-1) + aa*dcdt(i)*dTp(i);
    end
    Cond = Cond + ctm;
end
% [size(Cond),size(Temp),size(Pres)] % keyboard

Sdata.salin = sw_salt(Cond*10/sw_c3515, Temp, Pres); % Temp s/b IPTS-68?
Sdata.C_adj = Cond;
Sdata.T_adj = Temp;
