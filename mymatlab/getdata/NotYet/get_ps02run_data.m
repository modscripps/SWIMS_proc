function [SWd,ADd,EXd] = get_ps02run_data(section_num)
% Usage: [SWd,ADd,EXd] = get_ps02run_data(section_num);
% Gather SWIMS, ADCP data for a ps02 section, organized by run numbers.
%   section_num = 1 or 2 (1=HC2xN, 2=HC2xS);
% OUTPUT structures:
%   SWd(run) has SWIMS data (including along-thalweg/crossing distances) 
%   ADd(run) has ADCP data (including along-thalweg/crossing distances), and
%   EXd contains ancillary data for the section.
% Data are omitted for portions where reference distance is non-monotonic.
% DPW - Feb 2003
crz = 'ps02';
set_swims_paths

%% Edit to specify which section: (1=HC2xN, 2=HC2xS)
iSN = section_num;

%% Gather parameters for section from Run Index
load(fullfile(swimsindex,'SWIMS_ps02_Runs.mat'))
YD_BEG = Sec_params(iSN).start_yday; YD_END = Sec_params(iSN).end_yday;
TWfil = fullfile(swimsthalweg, Sec_params(iSN).Thalweg_file); 
% set up variable names for thalweg reference
TWstr = Sec_params(iSN).TWfil_structure;
TWvars = {[TWstr '.lon_hires'], [TWstr '.lat_hires'], ...
        [TWstr '.dist_hires'], [TWstr '.ang_hires']};
XLims = Sec_params(iSN).TWdist_lims; % round down/up 0.1 km
XLims = [floor(XLims(1)*10)/10, ceil(XLims(2)*10)/10];
TWang = Sec_params(iSN).Ref_values; % has fixed angle for crossings
if length(TWang)>2, TWang = TWang(3)-90; end % to point down-thalweg
ZLims = [0 max(Runs(iSN).maxdepth) + 5]; % deepest run + 5m (OK for ADCP)
ADfrac = 0.93; % WH ADCP to 93% of bottom
if iSN>0, ADfrac = 0.84; end % BB ADCP to 85% of bottom

%% Save some info
EXd.SECpre = Sec_params(iSN).Tag; % Prefix for output files
EXd.XLims = XLims;
EXd.ZLims = ZLims;
EXd.TWang = TWang;

%% Get all the section data (-/+ 14.4 min)
GP=get_gps_data(YD_BEG-.011, YD_END+.011, ...
    fullfile(swimsindex,'GPS_ps02_matfiles.mat'), fullfile(swimsmatdata,'GPS'));
SW=get_swims_data(YD_BEG-.01, YD_END+.01, ...
    fullfile(swimsindex,'SWIMS_ps02_gridfiles.mat'), swimsgridded); % , ...
%    {'starts','ends','th1','s1','sgth1','updown','etc...'});
AD=get_adcp_data(YD_BEG-.01, YD_END+.01, ...
    fullfile(swimsindex,'ADCP_ps02_matfiles.mat'), fullfile(swimsmatdata,'ADCP'), ...
    {'bottomBT','u_wat','v_wat','w_wat','ldayFP','ldayLP'});

%% Put in filler if no overturn stats available
if ~isfield(SW, 'eps1')
    SW.eps1 = NaN*SW.th1; SW.krho1 = NaN*SW.th1; SW.ovt1_Th = NaN*SW.th1;
end
%% Eliminate out-of-depth-range grid rows; put NaNs in contaminated ADCP bins
ix = find(SW.z_ctd < ZLims(1) | SW.z_ctd > ZLims(2));
if ~isempty(ix)
    SW.th1(ix,:) = []; SW.s1(ix,:) = []; SW.sgth1(ix,:) = []; SW.z_ctd(ix) = [];
    SW.dox(ix,:) = []; SW.flu(ix,:) = []; SW.obs(ix,:) = [];
    SW.eps1(ix,:) = []; SW.krho1(ix,:) = []; SW.ovt1_Th(ix,:) = [];
end
% Put in NaNs for w_wat if not saved in pre-042902 files
if ~isfield(AD,'w_wat')
    AD.w_wat = AD.u_wat*NaN;
end
ix = find(AD.z_adcp < ZLims(1) | AD.z_adcp > ZLims(2));
if ~isempty(ix)
    AD.u_wat(ix,:) = []; AD.v_wat(ix,:) = [];
    AD.w_wat(ix,:) = [];AD.z_adcp(ix) = [];
end
ix = find(isnan(AD.bottomBT)); % lost bottom-tracking
AD.u_wat(:,ix) = NaN; AD.v_wat(:,ix) = NaN; AD.w_wat(:,ix) = NaN;
% AD.u_wat(1,:) = NaN; AD.v_wat(1,:) = NaN; % first bin bogus (?)
for id=1:length(AD.bottomBT)
    ix = find(AD.z_adcp > AD.bottomBT(id)*ADfrac);
    if ~isempty(ix)
        AD.u_wat(ix,id) = NaN; AD.v_wat(ix,id) = NaN;
        AD.w_wat(ix,id) = NaN;
    end
end

%% Determine times,locations at start,middle,end of SWIMS/ADCP profiles
samps = SW.ends - SW.starts; % # of 24-Hz measurements during profile
prof_time = (samps/24) / 86400; % Duration (in days)
SW.yd_stop = SW.yday + prof_time; % end of profile
SW.yd_mid = (SW.yday + SW.yd_stop) / 2; % mid point of profile
% Midnight problem for ADCP GPS-days, quick fix (other ydays fixed in get-routines)
ix=find(AD.ldayFP<AD.ldayFP(1)); AD.ldayFP(ix)=AD.ldayFP(ix)+1;
ix=find(AD.ldayLP<AD.ldayLP(1)); AD.ldayLP(ix)=AD.ldayLP(ix)+1;
AD.ldayMP = (AD.ldayFP+AD.ldayLP)/2; % mid-profile yearday
% Interpolate positions from 1-sec GPS data
SW.lon_mid = interp1(GP.sattime,GP.lon,SW.yd_mid);
SW.lat_mid = interp1(GP.sattime,GP.lat,SW.yd_mid);
SW.lon_start = interp1(GP.sattime,GP.lon,SW.yday);
SW.lat_start = interp1(GP.sattime,GP.lat,SW.yday);
SW.lon_stop = interp1(GP.sattime,GP.lon,SW.yd_stop);
SW.lat_stop = interp1(GP.sattime,GP.lat,SW.yd_stop);
%
AD.lonMP = interp1(GP.sattime,GP.lon,AD.ldayMP);
AD.latMP = interp1(GP.sattime,GP.lat,AD.ldayMP);

clear GP

%% Now, loop thru Runs on the Section, gathering data into STRuctures(iRN).
for iRN=1:length(Runs(iSN).yd_start)
    % Get indices for profiles in time range
    YD_b = Runs(iSN).yd_start(iRN); YD_e = Runs(iSN).yd_end(iRN);
    iSyd = find(SW.yd_mid>YD_b & SW.yd_mid<YD_e); % exclusive limits
    iAyd = find(AD.ldayMP>YD_b & AD.ldayMP<YD_e); % may have <1min while turning
    iBTyd = [iAyd(1)-1 iAyd iAyd(end)+1]; % extra on ends for plotting bottom
    if iBTyd(1)<1, iBTyd(1)=[]; end
    if iBTyd(end)>length(AD.ldayMP), iBTyd(end)=[]; end
        
    % Compute thalweg references
    RFb = thalweg_refer(SW.lon_start(iSyd), SW.lat_start(iSyd), ...
        TWfil, TWvars, XLims(1)-.2, XLims(2)+.2);
    RFm = thalweg_refer(SW.lon_mid(iSyd), SW.lat_mid(iSyd), ...
        TWfil, TWvars, XLims(1)-.2, XLims(2)+.2);
    RFe = thalweg_refer(SW.lon_stop(iSyd), SW.lat_stop(iSyd), ...
        TWfil, TWvars, XLims(1)-.2, XLims(2)+.2);
    RFad = thalweg_refer(AD.lonMP(iAyd), AD.latMP(iAyd), ...
        TWfil, TWvars, XLims(1)-.2, XLims(2)+.2);
    % Substitute constant TW direction, if specified for Section
    if ~isnan(TWang), RFad.thal_ang(:) = TWang; end
    RFbt = thalweg_refer(AD.lonMP(iBTyd), AD.latMP(iBTyd), ...
        TWfil, TWvars, XLims(1)-.2, XLims(2)+.2);
    
    TWsign = 1; % positive TW direction
    if (RFe.thal_dist(end) - RFe.thal_dist(1)) < 0
        TWsign = -1; % negative TW direction
    end
    % Exclude profiles where boat's thalweg-relative direction changed
    xL = RFm.thal_dist(1) - TWsign; % distance at last good midpoint
    iSok = []; dz =  diff(SW.z_ctd(1:2)); % SWIMS; also check profile thickness
    YD_rev = []; % record periods of direction reversal
    for ip=1:length(iSyd)
        norev = ( ( RFm.thal_dist(ip)-RFb.thal_dist(ip) )*TWsign > 0 ...
                & ( RFe.thal_dist(ip)-RFm.thal_dist(ip) )*TWsign > 0 ...
                & ( RFm.thal_dist(ip)-xL )*TWsign > 0 );
        if norev & length( find(~isnan(SW.th1(:,iSyd(ip)))) ) * dz > 7
            iSok = [iSok ip];
            xL = RFm.thal_dist(ip);
        end
        if ~norev
            YD_rev = [YD_rev, [SW.yday(iSyd(ip)); SW.yd_stop(iSyd(ip))] ];
        end
    end
    iSyd = iSyd(iSok); % indices into SWIMS grids,ydays (iSok into RF?.var(s))
    % ADCP
    xL = RFad.thal_dist(1) - TWsign;
    iAok = [];
    for ip=1:length(iAyd)
        if ( RFad.thal_dist(ip)-xL )*TWsign > 0
            iAok = [iAok ip];
            xL = RFad.thal_dist(ip);
        end
    end
    iBok = iAok + 1; % Add points past ends for BT bottom depth, if dir. okay
    if iAyd(1)>iBTyd(1) & ...
            ( RFad.thal_dist(1)-RFbt.thal_dist(1) )*TWsign > 0
        iBok = [1 iBok];
    end
    if iAyd(end)<iBTyd(end) & ...
            ( RFbt.thal_dist(end)-RFad.thal_dist(iAok(end)) )*TWsign > 0
        iBok = [iBok length(RFbt.thal_dist)];
    end
    iAyd = iAyd(iAok); % indices into ADCP grids,ydays (iAok into RFad.var(s))
    iBTyd = iBTyd(iBok); % indices into ADCP BT,ydays (iBok into RFbt.var(s))
    
    %% Gather Run data
    SWd(iRN).run=iRN; SWd(iRN).YD_b=YD_b; SWd(iRN).YD_e=YD_e;
    SWd(iRN).YD_rev=YD_rev; % excluded reversal periods
    % yeardays, TW distances at beginning, middle, end of each 'profile'
    SWd(iRN).yday_bme = [SW.yday(iSyd); SW.yd_mid(iSyd); SW.yd_stop(iSyd)];
    SWd(iRN).TWdist_bme = ...
        [RFb.thal_dist(iSok); RFm.thal_dist(iSok); RFe.thal_dist(iSok)];
    SWd(iRN).z_ctd = SW.z_ctd;
    SWd(iRN).updown = SW.updown(iSyd);
    SWd(iRN).th1 = SW.th1(:,iSyd);
    SWd(iRN).s1 = SW.s1(:,iSyd);
    SWd(iRN).sgth1 = SW.sgth1(:,iSyd);
    SWd(iRN).dox = SW.dox(:,iSyd);
    SWd(iRN).flu = SW.flu(:,iSyd);
    SWd(iRN).obs = SW.obs(:,iSyd);
    SWd(iRN).eps1 = SW.eps1(:,iSyd);
    SWd(iRN).krho1 = SW.krho1(:,iSyd);
    SWd(iRN).ovt1_Th = SW.ovt1_Th(:,iSyd);
    
    %% Transform ADCP into along,across TW (acw from along)
    ADd(iRN).yday_MP = AD.ldayMP(iAyd);
    ADd(iRN).TWdist_MP = RFad.thal_dist(iAok);
    ADd(iRN).TWang = RFad.thal_ang(iAok);
    ADd(iRN).z_adcp = AD.z_adcp;
    ADd(iRN).u_wat = AD.u_wat(:,iAyd);
    ADd(iRN).v_wat = AD.v_wat(:,iAyd);
    ADd(iRN).w_wat = AD.w_wat(:,iAyd);
    ca = ones(size(AD.z_adcp)) * cos(ADd(iRN).TWang * pi/180);
    sa = ones(size(AD.z_adcp)) * sin(ADd(iRN).TWang * pi/180);
    ADd(iRN).u_along = ca.*AD.u_wat(:,iAyd) + sa.*AD.v_wat(:,iAyd);
    ADd(iRN).u_cross = -sa.*AD.u_wat(:,iAyd) + ca.*AD.v_wat(:,iAyd);
    % Bottom trace
    ADd(iRN).dist_BT = RFbt.thal_dist(iBok);
    ADd(iRN).bot_BT = AD.bottomBT(iBTyd);
    
end
