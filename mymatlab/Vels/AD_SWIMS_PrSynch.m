% AD_SWIMS_PrSynch.m
% cruise='ML04'; %cruise = 'bs03';
% SWIMSfoldersPC
% cruise, yd_b, yd_e should already be specified before calling this script

AD_SWIMS_PrSynch_params

SynchUPDN = 1; % assume UP and DN are always synch'd
iNS=find( isnan(EnsDNtoUP(:,2)) );
if ~isempty(iNS) % check if profile is in an asynch'd period
    ydl = [ EnsDNtoUP(:,1); inf ]; % yearday limits (CTD time reference)
    for i=1:length(iNS)
        if ( yd_b>=ydl(iNS(i)) || yd_e>=ydl(iNS(i)) ) && ...
            ( yd_b<=ydl(iNS(i)+1) || yd_e<=ydl(iNS(i)+1) )
            SynchUPDN = 0; % some portion of the profile IS un-synch'd
            break;
        end
    end
end
clear iNS ydl

ex_yday = 1.5/1440; % extra data at ends (enough to cover synch offsets, at least)

crz=cruise;
set_swims_paths

SP=get_SWIMS_SciData(yd_b-ex_yday, yd_e+ex_yday, ...
    fullfile(swimsindex,['CTD_' cruise '_matfiles.mat']), ...
    fullfile(swimsmatdata,'CTD'), 1);
%LD = get_WinchLine_data(yd_b-ex_yday, yd_e+ex_yday, ...
LD = get_WinchLineSIO_data(yd_b-ex_yday, yd_e+ex_yday, ...
    fullfile(swimsindex,['LD_' cruise '_matfiles.mat']), ...
    fullfile(swimsmatdata,'LD') );

clear SW
SW.Pr = medfilt1(SP.Pr,25);
SW.yday = SP.yday_adj;
SW.Pitch = SP.Pitch;
SW.Roll = SP.Roll;
switch lower(cruise)
    case {'philex08','wa_nliw_apr2013'} % Short period when T1,C1 were clogged, so ...
        SW.Temp = SP.T2;
        SW.Sal = SP.S2;
        SW.SalAlt = SP.S1;
    otherwise
        SW.Temp = SP.T1;
        SW.Sal = SP.S1;
        SW.SalAlt = SP.S2;
end
% 'fix' obvious salinity errors
ix = find(isnan(SW.Sal) | SW.Sal>40 | SW.Sal<0); % was 0.040 CU, now 40 PSU
SW.Sal(ix) = SW.SalAlt(ix);
if ~isreal(SW.Sal)
        ii = find(imag(SW.Sal));
        SW.Sal = real(SW.Sal);
        SW.Sal(ii) = NaN;
end
ix = find(isnan(SW.Sal) | SW.Sal>40 | SW.Sal<0);
iy = find(~isnan(SW.Sal) & SW.Sal<=40 & SW.Sal>=0);
sa = 35; % 0.035 CU
if ~isempty(iy)
    sa = mean(SW.Sal(iy));
end
SW.Sal(ix) = sa; % okay for now   

clear SP
SW.Svel = medfilt1(sw_svel(SW.Sal,SW.Temp,SW.Pr),25);
ptsC = 24; % number of 24-Hz points for center difn fall rate (here = 1 s)
%SW.Wpr = 100*(SW.Pr(ptsC+1:end)-SW.Pr(1:end-ptsC)) / (ptsC/24); % MPa
SW.Wpr = 1*(SW.Pr(ptsC+1:end)-SW.Pr(1:end-ptsC)) / (ptsC/24); % dbar
SW.Wdy = (SW.yday(ptsC+1:end)+SW.yday(1:end-ptsC)) / 2;

clear VelDN VelUP % get pre-processed Velocity data, in leveled-ADCP coords
VelDN = get_ADupdn_data(yd_b-ex_yday, yd_e+ex_yday, ...
    fullfile(swimsindex,['VelDN_' cruise '_matfiles.mat']), ...
    fullfile(swimsmatdata,'VelDN') );
VelUP = get_ADupdn_data(yd_b-ex_yday, yd_e+ex_yday, ...
    fullfile(swimsindex,['VelUP_' cruise '_matfiles.mat']), ...
    fullfile(swimsmatdata,'VelUP') );

%VelUP = []; % for wa_nliw_apr2013 (and maybe others), ADUP pitch/roll jumps
clear xx; % assign empty structure fields if either VelUP or VelDN absent
if isempty(VelUP)
    if isempty(VelDN)
        disp('No VelUP and VelDN data, exit.')
        return
    end
    xx = VelDN;
    xx(2).yday = NaN; xx(2).ens_no = NaN; xx(2).SWIMS_headM = NaN;
    clear VelUP
    VelUP = xx(2); % assigns same fields as VelDN, all empty
end
if isempty(VelDN)
    xx = VelUP;
    xx(2).yday = NaN; xx(2).ens_no = NaN; xx(2).SWIMS_headM = NaN;
    clear VelDN
    VelDN = xx(2); % assigns same fields as VelUP, all empty
    nBt = {'uBT','vBT','wBT','werrBT','rangeBT','rangeBTmin'};
    for i=1:length(nBt)
        VelDN.(nBt{i}) = NaN;
    end
end
clear xx

% load Vel_Proc.mat
% ADDN = VelDN;
% ADUP = VelUP;
% clear VelDN VelUP

% Make adjustment to true heading here
VelDN.headT = mod(VelDN.SWIMS_headM + MagDcl, 360);
VelUP.headT = mod(VelUP.SWIMS_headM + MagDcl, 360);

% Synch UP ensemble numbers to DN ones
VelUP.ens_orig = VelUP.ens_no; % save for reference
if SynchUPDN
    ix = find( EnsDNtoUP(:,2) ); % non-zero offsets
    ydl = [ EnsDNtoUP(:,1); inf ]; % yearday limits 
    for i=1:length(ix)
        ie = find( VelUP.yday>=ydl(ix(i)) & VelUP.yday<ydl(ix(i)+1) );
        VelUP.ens_no(ie) = VelUP.ens_no(ie) - EnsDNtoUP(ix(i),2);
    end
end
clear mb ix
% Here, make sure ensemble no.s are monotonic 
%% (This shouldn't really happen during a depth cycle, unless the
%%  CTD deck box was turned off and on to reset something)
xeD = [find(diff(VelDN.ens_no)<1) length(VelDN.ens_no)];
xeU = [find(diff(VelUP.ens_no)<1) length(VelUP.ens_no)];
mxeD = NaN; mxeU = NaN;
for i=1:length(xeD)-1
    mxeD(i) = VelDN.ens_no(xeD(i)); % keep these to use for VelUP to keep ens_nos synched
    % force to be monotonic
    VelDN.ens_no([xeD(i)+1:xeD(i+1)]) = VelDN.ens_no([xeD(i)+1:xeD(i+1)]) + mxeD(i)+50;
end
if SynchUPDN % UP synched to DN
    for i=1:min(length(xeU),length(xeD))-1
        % force to be monotonic (use VelDN 'offsets' to keep synched)
        VelUP.ens_no([xeU(i)+1:xeU(i+1)]) = VelUP.ens_no([xeU(i)+1:xeU(i+1)]) + mxeD(i)+50;
    end
    if length(xeU) > length(xeD)
        VelUP.ens_no([xeU(length(xeD))+1:end]) = inf; % beyond last VelDN ensemble
    end
else % Non-synched, do UP by itself
    for i=1:length(xeU)-1
        mxeU(i) = VelUP.ens_no(xeU(i)); % really don't need these, but why not
        % force to be monotonic
        VelUP.ens_no([xeU(i)+1:xeU(i+1)]) = VelUP.ens_no([xeU(i)+1:xeU(i+1)]) + mxeU(i)+50;
    end
end
   
% Ranges of DN, UP ensemble numbers, initialize
en_bD = NaN; en_eD = NaN;
en_bU = NaN; en_eU = NaN;

% Adjust VelDN (water ping) time to agree with SWIMS (CTD)
VelDN.yday_DN = VelDN.yday + ...
    interp1(DelDNtoSW(:,1), DelDNtoSW(:,2), VelDN.yday, 'linear','extrap');
% iDN specifies indices for VelDN pings in yearday (and ensemble) range
iDN = find(VelDN.yday_DN>=yd_b & VelDN.yday_DN<=yd_e);
if ~isempty(iDN)
    en_bD = min(VelDN.ens_no(iDN));
    en_eD = max(VelDN.ens_no(iDN));
    % Prepare yearday, pressure, soundvel, synched to CTD; also CTD-fallrate
    ix = find(~isnan(SW.Svel));
    VelDN.pr_DN = interp1(SW.yday, SW.Pr, VelDN.yday_DN);
    VelDN.svel_DN = interp1(SW.yday(ix), SW.Svel(ix), VelDN.yday_DN);
    VelDN.w_DN = interp1(SW.Wdy, SW.Wpr, VelDN.yday_DN);
    % BT ping is approx 0.4-s before VelDN water ping
    VelDN.yday_BT = VelDN.yday_DN + DelDNtoBT; 
    VelDN.pr_BT = interp1(SW.yday, SW.Pr, VelDN.yday_BT);
    VelDN.svel_BT = interp1(SW.yday(ix), SW.Svel(ix), VelDN.yday_BT);
    %VelDN.BOT = VelDN.pr_BT*100 + VelDN.rangeBT; % first guess, no soundspeed adjustment
    VelDN.BOT = VelDN.pr_BT + VelDN.rangeBT; % pr/dbar,     "           "
end

% Adjust VelUP (water ping) time to agree with SWIMS (CTD)
VelUP.yday_UP = NaN * VelUP.ens_no; % in case some ensembles (pings) are missing
if SynchUPDN && ~isempty(iDN) % synch VelUP using ensemble numbers
    iUP = find(VelUP.ens_no>=en_bD & VelUP.ens_no<=en_eD); 
    VelUP.yday_UP(iUP) = ...
        interp1(VelDN.ens_no(iDN), VelDN.yday_DN(iDN), VelUP.ens_no(iUP)) + DelDNtoUP;
else % Need array DelUPtoSW(:,1:2) to map from ADUP to CTD time
    VelUP.yday_UP = VelUP.yday + ...
        interp1(DelUPtoSW(:,1), DelUPtoSW(:,2), VelUP.yday, 'linear','extrap');
    iUP = find(VelUP.yday_UP>=yd_b & VelUP.yday_UP<=yd_e);
end
% iUP specifies indices for VelUP pings in yearday (and ensemble) range
if ~isempty(iUP)
    en_bU = min(VelUP.ens_no(iUP));
    en_eU = max(VelUP.ens_no(iUP));
    % Prepare yearday, pressure, soundvel, synched to CTD; also CTD-fallrate
    ix = find(~isnan(SW.Svel));
    VelUP.pr_UP = interp1(SW.yday, SW.Pr, VelUP.yday_UP);
    VelUP.svel_UP = interp1(SW.yday(ix), SW.Svel(ix), VelUP.yday_UP);
    VelUP.w_UP = interp1(SW.Wdy, SW.Wpr, VelUP.yday_UP);
end

if ~SynchUPDN
    return % just to give debug point
end
return
