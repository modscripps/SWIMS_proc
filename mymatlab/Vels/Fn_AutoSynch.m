% AD_AutoSynch.m - For a specified time period, estimate ADCP ensemble
%   offset (ADUP-ADDN) and clock offset (CTD-ADDN) for SWIMS data
function [DtoU, DNtoCTD] = Fn_AutoSynch(yb, ye)

crz='ArcticMix';
swimsfolders
SG=get_swims_data(yb,ye, ...
    fullfile(swimsindex, 'SWIMS_ArcticMix_gridfiles.mat'), swimsgridded);
if isempty(SG)
    DtoU = [];
    DNtoCTD = [];
    return
end
AU = get_ADupdn_data(yb-1/1440,ye+1/1440, ...
    fullfile(swimsindex,'ADUP_ArcticMix_matfiles.mat'),fullfile(swimsmatdata,'ADUP'));

AD = get_ADupdn_data(yb-1/1440,ye+1/1440, ...
    fullfile(swimsindex,'ADDN_ArcticMix_matfiles.mat'),fullfile(swimsmatdata,'ADDN'));

load(fullfile(swimsmatdata,'ADUP-EnsYdays.mat'))
ix = find(EnsYdays.PCyday>=yb & EnsYdays.PCyday<=ye);

AUeyd.ens_no = EnsYdays.ens_no(ix);
AUeyd.RDyday = EnsYdays.RDyday(ix);
AUeyd.PCyday = EnsYdays.PCyday(ix);
clear EnsYdays
% attempt to remove errant outliers
ensMD = median(AUeyd.ens_no);
RmPyd = (AUeyd.RDyday-AUeyd.PCyday)*86400;
rpMD = median(RmPyd);

if isempty(AUeyd.PCyday)
    DtoU = [];
    DNtoCTD = [];
    return
end
% unlikely ensemble numbers and ADUP-CTD time difn's:
mxct = (AUeyd.PCyday(end)-AUeyd.PCyday(1))*86400 /1; %max at 1-s/ens
ix = find(AUeyd.ens_no<ensMD-mxct | AUeyd.ens_no>ensMD+mxct | ...
    abs(RmPyd-rpMD)>60);
AUeyd.ens_no(ix)=[];
AUeyd.RDyday(ix)=[];
AUeyd.PCyday(ix)=[];
RmPyd(ix)=[];

load(fullfile(swimsmatdata,'ADDN-EnsYdays.mat'))
ix = find(EnsYdays.PCyday>=yb & EnsYdays.PCyday<=ye);
ADeyd.PCyday = EnsYdays.PCyday(ix);
ADeyd.RDyday = EnsYdays.RDyday(ix);
ADeyd.ens_no = EnsYdays.ens_no(ix);
clear EnsYdays
ensMD = median(ADeyd.ens_no);
RmPyd = (ADeyd.RDyday-ADeyd.PCyday)*86400;
rpMD = median(RmPyd);
% unlikely ensemble numbers and ADDN-CTD time difn's:
mxct = (ADeyd.PCyday(end)-ADeyd.PCyday(1))*86400 /1; %max at 1-s/ens
ix = find(ADeyd.ens_no<ensMD-mxct | ADeyd.ens_no>ensMD+mxct | ...
    abs(RmPyd-rpMD)>60);
ADeyd.ens_no(ix)=[];
ADeyd.RDyday(ix)=[];
ADeyd.PCyday(ix)=[];
RmPyd(ix)=[];

% Now, estimate offsets. Verify graphically using Run_SWIMS_sync.m (turned
% off plotting for looping through)
%   For now, enter as EnsDNtoUP and DelDNtoSW in AD_SWIMS_PrSynch_params.m
ensUonD = interp1(AUeyd.PCyday,AUeyd.ens_no,ADeyd.PCyday);
%figure;
%plot(ADeyd.PCyday, ensUonD-ADeyd.ens_no, '.r-'),grid on
ix = find(~isnan(ensUonD-ADeyd.ens_no));
% DtoU should be within +/- 2 ens, minimize ADUP.pitch-ADDN.roll to verify:
DtoU = round(median(ensUonD(ix)-ADeyd.ens_no(ix))); % EnsDNtoUP
%figure;
%plot(ADeyd.PCyday,(ADeyd.RDyday-ADeyd.PCyday)*86400,'.r-'),grid on
% DNtoCTD should be close to offset, but seems ~1.7-s to long
%   according to graphic comparison of fall/rise rates via AD_CTD_Synch.m 
DNtoCTD = -nanmax((ADeyd.RDyday-ADeyd.PCyday)*86400); % DelDNtoSW
end