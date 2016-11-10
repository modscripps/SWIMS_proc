function Vel = Get_ADCP_SwimsUP_info(RDIname, filetype, savetype);
% Vel = Get_ADCP_SwimsUP_info(RDIname, filetype, savetype);
%	Read in RDI file(s), put requested results in structure = Vel
%  filetype: 0=process single file=RDIname, 1=multi LTA files
%		starting with RDIname (more later);
%	savetype: 0=just depth,time,position,orientation vectors,
%		(No profile data for 'info' version)
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
%MagDcl = 19.5 ; % for Puget Sound, Hood Canal, approx
MagDcl = 10.5 ; % for Oahu, approx

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
   raw2mat_info_V3(RDIname,TMPfil, [], {'.mat'});
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
% transformation matrix, beam-to-adcp
TXF=  [1.4619   -1.4619    0.0000    0.0000;
  0.0000    0.0000   -1.4619    1.4619;
  0.2661    0.2661    0.2661    0.2661;
  1.0337    1.0337   -1.0337   -1.0337];

UpDown = 1; % -1=down, +1=upward

Vel.ens_no = ensemble_number + ensMSB*65536;

% Some will be the same for all ensembles:
Vel.z_adcp =  - UpDown*(dis1/100 + [0:nbins-1]*binlen/100); % meters
Vel.p_adcp = Vel.z_adcp/100; % MPa
Vel.pulselen = pulselen;

% Others will vary:
Vel.yday = yearday(day,month,year,hour,minute,second+hundreths/100);

iSV=1:length(Vel.yday);
%iSV=find(Vel.yday>101.7450 & Vel.yday<101.7510); % for PS apr 22,2002 dn/up test
Vel.yday = Vel.yday(iSV);
Vel.ens_no = Vel.ens_no(iSV);

% Orientation data (Internal for SWIMS)
Vel.heading = heading(iSV)/100;
Vel.pitch = pitch(iSV)/100;
Vel.roll = roll(iSV)/100;
% Find angles in SWIMS frame: heading toward nose, pitch nose up, roll to right;
% Assumes forward beams are 3^1 for upward, 1^3 for downward ADCP
Vel.SWIMS_headT = mod(Vel.heading + UpDown*45 + MagDcl, 360);
Adir = atan2(Vel.pitch,Vel.roll) - UpDown*pi/4;
Vel.SWIMS_pitch = sqrt(Vel.pitch.^2+Vel.roll.^2).*sin(Adir);
Vel.SWIMS_roll = sqrt(Vel.pitch.^2+Vel.roll.^2).*cos(Adir);

%% BT - range and velocity (downward only)
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
    
    %% Compute East,North BT vels (of SWIMS)
    BTdn = TXF*Vel.btvel_bm;
    ANGdn = 90-Vel.SWIMS_headT;
    Vel.VbtE = cos((ANGdn-45)*pi/180).*(-BTdn(1,:)) + ...
        sin((ANGdn-45)*pi/180).*(-BTdn(2,:));
    Vel.VbtN = -sin((ANGdn-45)*pi/180).*(-BTdn(1,:)) + ...
        cos((ANGdn-45)*pi/180).*(-BTdn(2,:));
    Vel.VbtZ = -BTdn(3,:);
    
end % UpDown<0 for BT

return



