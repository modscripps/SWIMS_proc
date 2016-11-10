function Vel = Get_ADCP_prof(RDIname, filetype, savetype);
% Vel = Get_ADCP_prof(RDIname, filetype, savetype);
%	Read in RDI file(s), put requested results in structure = Vel
%  filetype: 0=process single file=RDIname, 1=multi LTA files
%		starting with RDIname (more later);
%	savetype: 0=just u,v,w, depth,time,position vectors,
%		1=more stuff (%good, etc), 2=everything
%  Vel = structure with useful RDI data
% Dave W, 03-Apr-2001

if nargin<2
   filetype = 0; % single profile file (EG via serial Ensemble output)
end
if nargin<3
   savetype = 0;
end

%savetype=0; % for now

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
elseif filetype==4 % Transect Raw ADCP (multi-out later)
   raw2mat_V3(RDIname,TMPfil, {3}, {'.mat'});
end

% load in Matlab arrays, convert to useful units, save in structure
load(TMPfil)
delete([TMPfil '.mat']);
% Check if no valid data were found, exit if so
if isempty(ensemble_number)
    Vel = [];
    return
end

BADVEL = -32768;

Vel.ens_no = ensemble_number;

% Some will be the same for all ensembles:
Vel.z_adcp = dis1/100 + [0:nbins-1]*binlen/100; % meters
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
btrange = [btrange1; btrange2; btrange3; btrange4];
for i=1:4
   ib = find(btrange(i,:) == 0);
   if ~isempty(ib)
	   btrange(i,ib) = NaN; 
	end
end
Vel.bottomBT = nanmean(btrange)/100 + xducerdepth/10;

%
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
u_bot = btvel1/1000; v_bot = btvel2/1000; % m/s
w_bot = btvel3/1000; err_bot = btvel4/1000;

% Replace BT vels if missing or out-of-range with NAV, flag as non-BT vels
Vel.goodBT = ones(size(u_bot));
ixbt = find(isnan(u_bot+v_bot)); % used NAV for these
Vel.goodBT(ixbt) = 0;
BAM32 = (180/2^31); BAM16 = (180/2^15);% Nav speed,direction made good
Vel.lSMG = navSMG/1000; % speed made good
Vel.lDMG = navDMG * BAM16; % direction made good
if ~isempty(ixbt) % interpret in terms of bottom speed, to mimic BT sense
    u_bot(ixbt) = -Vel.lSMG(ixbt) .* cos( (90-Vel.lDMG(ixbt)) * pi/180 );
    v_bot(ixbt) = -Vel.lSMG(ixbt) .* sin( (90-Vel.lDMG(ixbt)) * pi/180 );
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
dayFP = yearday(navutcday,navutcmonth,navutcyear,0,0,navutctime1*1e-4);
dayLP = yearday(navutcday,navutcmonth,navutcyear,0,0,navutctime2*1e-4);
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

% Orientation data (external or internal)
Vel.heading = heading/100;
Vel.pitch = pitch/100;
Vel.roll = roll/100;

if savetype>0
   % Save more ...
   Vel.heading = heading;
end

return



