function Vel = Get_ADCP_profRaw(RDIname, filetype, savetype);
% Vel = Get_ADCP_profRaw(RDIname, filetype, savetype);
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

filetype = -1; % for Raw: single file with multiple ensembles

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

% Set up temporary file for matlab ADCP data
fiTMP = -1;
if 0==1
	for ix=0:1000
		TMPfil = ['TMP' num2str(ix)];
		fid = fopen([TMPfil '.mat'], 'r');
		if fid<0 % file doesn't exist yet
			fiTMP = fopen([TMPfil '.mat'], 'w');
			break
		else
			fclose(fid);
		end
	end
	if fiTMP<0
		error('Can not find available temporary file, exiting')
		return
	end
	fclose(fiTMP);
end

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

BADVEL = -32768;

UpDown = 1; % -1=down, +1=upward

xducerdepth = 47; % for PS apr 22,2002 dn/up test, holding at 47 m
Vel.ens_no = ensemble_number + ensMSB*65536;;

% Some will be the same for all ensembles:
Vel.z_adcp = xducerdepth - UpDown*(dis1/100 + [0:nbins-1]*binlen/100); % meters
Vel.p_adcp = Vel.z_adcp/100; % MPa
Vel.pulselen = pulselen;

% Others will vary:
Vel.yday = yearday(day,month,year,hour,minute,second+hundreths/100);

iSV=1:length(Vel.yday);
%iSV=find(Vel.yday>101.7450 & Vel.yday<101.7510); % for PS apr 22,2002 dn/up test
Vel.yday = Vel.yday(iSV);
Vel.ens_no = Vel.ens_no(iSV);

%% BT - range and velocity
if UpDown<0
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

end % UpDown<0 for BT

v1=v1(:,iSV);v2=v2(:,iSV);v3=v3(:,iSV);v4=v4(:,iSV);

%% Measured beam components
bad = find(v1==BADVEL); if bad, v1(bad) = NaN; end
bad = find(v2==BADVEL); if bad, v2(bad) = NaN; end
bad = find(v3==BADVEL); if bad, v3(bad) = NaN; end
bad = find(v4==BADVEL); if bad, v4(bad) = NaN; end
Vel.v1_bm = v1/1000; Vel.v2_bm = v2/1000; % m/s
Vel.v3_bm = v3/1000; Vel.v4_bm = v4/1000; % m/s

% Orientation data (Internal for SWIMS)
Vel.heading = heading(iSV)/100;
Vel.pitch = pitch(iSV)/100;
Vel.roll = roll(iSV)/100;

if savetype>0
   % Save more ...
   Vel.ec1_bm = e1; Vel.ec2_bm = e2; Vel.ec3_bm = e3; Vel.ec4_bm = e4;
   Vel.cor1_bm = cor1; Vel.cor2_bm = cor2; Vel.cor3_bm = cor3; Vel.cor4_bm = cor4;
end

return



