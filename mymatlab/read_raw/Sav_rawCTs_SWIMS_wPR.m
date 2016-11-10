function Raw_wPR = Sav_rawCTs_SWIMS_wPR(RawStr)
% usage: Raw_wPR = Sav_rawCTs_SWIMS_wPR(RawStr);
%   RawStr = Raw CTD data structure, as returned from Read_rawCTs_SWIMS.m
%       or Read_rawCTsHEX_SWIMS.m
%   RwPR = Structure of CTD raw data saved in CTD-*.mat files;
%       same as input structure, but with added fields of
%       calibrated pressure [dbar] and with SWIMStime converted from
%       datenum to yearday (for backwards compatability with pre-2015 data)
%   Dave W - 6/2015

% copy into output structure
Raw_wPR = RawStr; clear RawStr

year0 = datevec(Raw_wPR.SWIMStime(1)); year0 = year0(1);
yd0 = datenum(year0, 1, 1, 0, 0, 0);
Raw_wPR.SWIMStime = Raw_wPR.SWIMStime - yd0; % convert to yeardays

% compute and save pressure / dbar;
%   smooth temp.comp. for pressure over 10 sec window
smL=min(floor(length(Raw_wPR.tprCT)/2)-2,2401);
if ~mod(smL,2), smL = smL-1; end % must be odd
smC=ones(1,smL)/smL; smH=floor(smL/2);
tprSM = conv(Raw_wPR.tprCT,smC);
tprSM(1:smH)=[]; tprSM(1:smH)=tprSM(smH+1);
tprSM(end-smH+1:end)=[]; tprSM(end-smH+1:end)=tprSM(end-smH);
[prn, xx] = GetSWIMSConfig('pr', year0, Raw_wPR.SWIMStime(1));
[Pcoefs] = read_Pcal2_sbe_SWIMS(prn, year0, Raw_wPR.SWIMStime(1));
[TPRcoefs] = read_TPRcal_sbe_SWIMS(prn, year0, Raw_wPR.SWIMStime(1));
% compute temp, then pressure
Tpr = TPr_Swims(tprSM, TPRcoefs.m, TPRcoefs.b);
Raw_wPR.Pr = P_sbe2(Raw_wPR.pfreq/256, Tpr, Pcoefs);
clear Tpr tprSM
% %Fix up times here? - use modct,etc
%YDAYmin = Raw_wPR.SWIMStime(1); YDAYmax = Raw_wPR.SWIMStime(end);
