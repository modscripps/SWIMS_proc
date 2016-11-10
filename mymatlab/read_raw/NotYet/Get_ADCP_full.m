function Vel = Get_ADCP_full(RDIname, filetype, savetype);
% Vel = Get_ADCP_full(RDIname, filetype, savetype);
%	Read in RDI file(s), put requested results in structure = Vel
%  filetype: 0=process single file=RDIname, 1=multi LTA files
%		starting with RDIname (more later);
%	savetype: 0=just u,v,w, depth,time,position vectors,
%		1=more stuff (%good, etc), 2=everything
%  Vel = structure with useful RDI data
% Dave W, 03-Apr-2001
% Get_ADCP_profRaw.m modified for raw ping data received from SWIMS ADCPs
if nargin<2
   filetype = 0; % single profile file (EG via serial Ensemble output)
end
if nargin<3
   savetype = 0;
end

%filetype = -1; % for Raw: single file with multiple ensembles

Vel=[];

if filetype==0  % check that this is one complete file/profile
   fid=fopen(RDIname, 'r','ieee-le');
   if fid<0
      error(['Cannot open file = ' RDIname])
   end
   A = fread(fid,'char');
   frewind(fid);
   id = fread(fid,1,'uchar');
   src= fread(fid,1,'uchar');
   nbytes=fread(fid,1,'ushort');
   fclose(fid);
   if id~=127 | src~=127
      error([RDIname ' does not start with RDI header record']);
   end
   if length(A)-nbytes ~= 2
      error([RDIname ' is wrong size for single RDI profile']);
   end
end % of solo-file check

TMPfil='TMP0';

% Convert File(s)
if filetype < 1 % either: 0 = solo profile, or <0 = one file, multi-profs
   raw2mat_V3(RDIname,TMPfil, [], {'.mat'});
elseif filetype==1 % VmDas long-term avgs
   raw2mat_V3(RDIname,TMPfil, {6,'.LTA'}, {'.mat'});
elseif filetype==2 % VmDas short-term avgs
   raw2mat_V3(RDIname,TMPfil, {6,'.STA'}, {'.mat'});
elseif filetype==3 % VmDas Raw ADCP (multi-out later)
   raw2mat_V3(RDIname,TMPfil, {6,'.ENR'}, {'.mat'}); 
elseif filetype==4 % VmDas Raw ADCP (multi-out later)
   raw2mat_V3(RDIname,TMPfil, {6,'.ENX'}, {'.mat'}); 
elseif filetype==5 % Transect Raw ADCP (multi-out later)
   raw2mat_V3(RDIname,TMPfil, {3}, {'.mat'});
end

% load in Matlab arrays, convert to useful units, save in structure
load(TMPfil)
%keyboard
delete([TMPfil '.mat']);
% Check if no valid data were found, exit if so
if isempty(ensemble_number)
    Vel = [];
    return
end

BADVEL = -32768;

UpDown = -1; % -1=down, +1=upward

Vel.ens_no = ensemble_number + ensMSB*65536;

% Some will be the same for all ensembles:
Vel.z_adcp =  - UpDown*(dis1/100 + [0:nbins-1]*binlen/100); % meters
Vel.p_adcp = Vel.z_adcp/100; % MPa
Vel.pulselen = pulselen;

% Others will vary:
Vel.yday = yearday(day,month,year,hour,minute,second+hundreths/100);
%% BT - range and velocity
if isempty(btrange1) % If set for NOT bottom-tracking, put in NaNs (3-18-2003)
        ix = NaN*Vel.ens_no;
        btrange1 = ix; btrange2 = ix; btrange3 = ix; btrange4 = ix;
        btvel1 = ix; btvel2 = ix; btvel3 = ix; btvel4 = ix;
end

iSV=1:length(Vel.yday);
%iSV=find(Vel.yday>101.7450 & Vel.yday<101.7510); % for PS apr 22,2002 dn/up test
Vel.yday = Vel.yday(iSV);
Vel.ens_no = Vel.ens_no(iSV);

% Orientation data (Internal for SWIMS)
Vel.heading = heading(iSV)/100;
Vel.pitch = pitch(iSV)/100;
Vel.roll = roll(iSV)/100;
Vel.hdg_std = stdhed(iSV)/1;
Vel.pitch_std = stdpitch(iSV)/10;
Vel.roll_std = stdroll(iSV)/10;

btrange = [btrange1(iSV); btrange2(iSV); btrange3(iSV); btrange4(iSV)];
for i=1:4
    ib = find(btrange(i,:) == 0);
    if ~isempty(ib)
        btrange(i,ib) = NaN; 
    end
end
Vel.btrange = btrange;
Vel.bottomBT = nanmean(btrange)/100 + xducerdepth/10;

%
btvel1=btvel1(iSV);btvel2=btvel2(iSV);btvel3=btvel3(iSV);btvel4=btvel4(iSV);
ib = find(btvel1==BADVEL);
if ~isempty(ib)
    btvel1(ib) = NaN; 
end
ib = find(btvel2==BADVEL);
if ~isempty(ib)
    btvel2(ib) = NaN;
end
ib = find(btvel3==BADVEL);
if ~isempty(ib)
    btvel3(ib) = NaN;
end
Vel.btvel_bm = [btvel1; btvel2; btvel3; btvel4] / 1000; % m/s
        
v1=v1(:,iSV);v2=v2(:,iSV);v3=v3(:,iSV);v4=v4(:,iSV);

%% Measured beam components
bad = find(v1==BADVEL); if bad, v1(bad) = NaN; end
bad = find(v2==BADVEL); if bad, v2(bad) = NaN; end
bad = find(v3==BADVEL); if bad, v3(bad) = NaN; end
bad = find(v4==BADVEL); if bad, v4(bad) = NaN; end
Vel.v1 = v1/1000; Vel.v2 = v2/1000; % m/s
Vel.v3 = v3/1000; Vel.v4 = v4/1000; % m/s

if savetype>0 % need to include iSV if other than all (columns)
    % Save more ...
    Vel.ec1_bm = e1; Vel.ec2_bm = e2; Vel.ec3_bm = e3; Vel.ec4_bm = e4;
    Vel.cor1_bm = cor1; Vel.cor2_bm = cor2; Vel.cor3_bm = cor3; Vel.cor4_bm = cor4;
    Vel.pg1 = pg1; Vel.pg2 = pg2; Vel.pg3 = pg3; Vel.pg4 = pg4;
    if savetype>=2  & ~isempty(navlat1) % implement this later
        u_bot = btvel1/1000; v_bot = btvel2/1000; % m/s
        w_bot = btvel3/1000; err_bot = btvel4/1000;
        
        % Replace BT vels if missing or out-of-range with NAV, flag as non-BT vels
        Vel.goodBT = ones(size(u_bot));
        ixbt = find(isnan(u_bot+v_bot)); % used NAV for these
        Vel.goodBT(ixbt) = 0;
        BAM32 = (180/2^31); BAM16 = (180/2^15);% Nav speed,direction made good
        Vel.lSMG = navSMG/1000; % speed made good (bkwd dif = <E_n>-<E_n-1>)
        Vel.lDMG = navDMG * BAM16; % direction made good
        Vel.NavSpd = navavgspd/1000; % avg navg speed (in ensemble = ?)
        Vel.NavDirT = navavgtrackT * BAM16; % avg navg direction True
        Vel.NavDirM = navavgtrackM * BAM16; % avg navg direction Mag
        Vel.Nav_U = Vel.NavSpd .* cos( (90-Vel.NavDirT) * pi/180 );
        Vel.Nav_V = Vel.NavSpd .* sin( (90-Vel.NavDirT) * pi/180 );
        if ~isempty(ixbt) % interpret in terms of bottom speed, to mimic BT sense
            u_bot(ixbt) = -Vel.lSMG(ixbt) .* cos( (90-Vel.lDMG(ixbt)) * pi/180 );
            v_bot(ixbt) = -Vel.lSMG(ixbt) .* sin( (90-Vel.lDMG(ixbt)) * pi/180 );
            % Use avg within ensemble instead?  Then set:
            % u_bot(ixbt) = -Vel.Nav_U(ixbt);
            % v_bot(ixbt) = -Vel.Nav_V(ixbt); 
        end
        
        %% Measured and BT-corrected water velocity
        bad = find(v1==BADVEL); if bad, v1(bad) = NaN; end
        bad = find(v2==BADVEL); if bad, v2(bad) = NaN; end
        bad = find(v3==BADVEL); if bad, v3(bad) = NaN; end
        bad = find(v4==BADVEL); if bad, v4(bad) = NaN; end
        u_rel = v1/1000; v_rel = v2/1000; % m/s
        w_rel = v3/1000; err_rel = v4/1000; % m/s
        %
        [m,n] = size(u_rel);
        if ~isempty(u_bot) 
            uno = ones(m,1);Vel.u_wat = u_rel - uno*u_bot; Vel.v_wat = v_rel - uno*v_bot;
            Vel.w_wat = w_rel - uno*w_bot;
            % Vessel velocity
            Vel.u_shipBT = -u_bot; Vel.v_shipBT = -v_bot; Vel.w_shipBT = -w_bot;
        else % should change to (e.g.) vessel=mean(u_rel,v_rel columns), wat=demeaned(u_rel,v_rel)
            Vel.u_wat=NaN*u_rel;Vel.v_wat=NaN*u_rel;Vel.w_wat=NaN*u_rel;
            Vel.u_shipBT=NaN;Vel.v_shipBT=NaN;Vel.w_shipBT=NaN;
        end
        
        %% Navigation data: FP,LP = First,Last Position during avg period
        BAM32 = (180/2^31); BAM16 = (180/2^15);
        Vel.latFP = navlat1 * BAM32;
        Vel.lonFP = navlon1 * BAM32;
        Vel.latLP = navlat2 * BAM32;
        Vel.lonLP = navlon2 * BAM32;
        % Nav UTC times are (1e-4) seconds since midnight
        dayFP = navutcday*NaN;
        dayLP = navutcday*NaN;
        iok = find(navutcday>=1 & navutcmonth>=1 & navutcyear>=1);
        x = yearday(navutcday(iok),navutcmonth(iok),navutcyear(iok), ...
            0,0,navutctime1(iok)*1e-4);
        dayFP(iok) = x;
        x = yearday(navutcday(iok),navutcmonth(iok),navutcyear(iok), ...
            0,0,navutctime2(iok)*1e-4);
        dayLP(iok) = x;
        % Check for RDI bug, not switching date for first fixes after midnight
        ixX = find(diff(dayFP)<0);
        while ~isempty(ixX)
            dayFP(ixX(1)+1) = dayFP(ixX(1)+1)+1;
            pause(1)
            ixX = find(diff(dayFP)<0);
        end
        ixX = find(diff(dayLP)<0);
        while ~isempty(ixX)
            dayLP(ixX(1)+1) = dayLP(ixX(1)+1)+1;
            pause(1)
            ixX = find(diff(dayLP)<0);
        end
        Vel.ldayFP = dayFP;
        Vel.ldayLP = dayLP;
        Vel.lPCtimeoff = navpctimeoff/1000; % PC-UTC(gps) seconds
    end % of "if 0 (savetype>=2)"
end

return



