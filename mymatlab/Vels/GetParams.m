function [EnsDNtoUP, DelDNtoSW] = GetParams()
% GetParams.m Returns EnsDNtoUp and DelDNtoSW arrays for use in
% AD_SWIMS_PrSynch_params.m 
%   Tests each half hour of cruise

sday = 238; %start of cruise
eday = 268; %end of cruise
inc = 0.05;
EnsDNtoUP = [];
DelDNtoSW = [];
oldDNtoCTD = 0;

for yb = sday:.001:eday-inc
    ye = yb+inc;
    [DtoU, DNtoCTD] = Fn_AutoSynch(yb, ye);
    if abs(DNtoCTD-oldDNtoCTD) > 1
        EnsDNtoUP = [EnsDNtoUP; yb DtoU];
        DelDNtoSW = [DelDNtoSW; yb DNtoCTD];
        oldDNtoCTD = DNtoCTD;
    end
end

