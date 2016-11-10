function Vel = Get_ADCP_SwimsDN(RDIname, filetype, savetype, ByteRange);
% Vel = Get_ADCP_SwimsDN(RDIname, filetype, savetype);
%	Read in RDI file(s), put requested results in structure = Vel
%  filetype: 0=process single file=RDIname, 1=multi LTA files
%		starting with RDIname (more later);
%	savetype: 0=just u,v,w, depth,time,position vectors,
%		1=more stuff (%good, etc), 2=everything;
%   ByteRange: range of byte numbers to process in input file (Aug 2015)
%  Vel = structure with useful RDI data - added LastByte in Aug 2015
%   to report how far file was actually processed
% Dave W, 03-Apr-2001
% Get_ADCP_profRaw.m modified for raw ping data received from SWIMS ADCPs
if nargin<2
   filetype = 0; % single profile file (EG via serial Ensemble output)
end
if nargin<3
   savetype = 0;
end
if nargin<4 || isempty(ByteRange)
    ByteRange = [0,inf];
end

filetype = -1; % hard-coded for SWIMS Raw: single file, multiple ensembles
%MagDcl = 19.5 ; % for Puget Sound, Hood Canal, approx
MagDcl = 10.5 ; % for Oahu, approx
%MagDcl = 3.5; % for Black Sea, approx (was set at 10.5 until yday>76.4)

Vel=[];

[fid,message]=fopen(filename,'r');

if fseek(fid, ByteRange(1), 'bof') < 0  % start reading here
    warning(['starting byte ' num2str(ByteRange(1)) ' past end-of-file: ',filename])
    fclose(fid);
    return
end
fclose(fid);

% if filetype==0  % check that this is one complete file/profile
%    fid=fopen(RDIname, 'r','ieee-le');
%    if fid<0
%       error(['Cannot open file = ' RDIname])
%    end
%    A = fread(fid,'char');
%    frewind(fid);
%    id = fread(fid,1,'uchar');
%    src= fread(fid,1,'uchar');
%    nbytes=fread(fid,1,'ushort');
%    fclose(fid);
%    if id~=127 | src~=127
%       error([RDIname ' does not start with RDI header record']);
%    end
%    if length(A)-nbytes ~= 2
%       error([RDIname ' is wrong size for single RDI profile']);
%    end
% end % of solo-file check

TMPfil='TMP0';

% Convert File(s)
if filetype < 1 % either: 0 = solo profile, or <0 = one file, multi-profs
   raw2mat_V4(RDIname,TMPfil, [], {'.mat'}, ByteRange);
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
% transformation matrix, beam-to-adcp
TXF=  [1.4619   -1.4619    0.0000    0.0000;
  0.0000    0.0000   -1.4619    1.4619;
  0.2661    0.2661    0.2661    0.2661;
  1.0337    1.0337   -1.0337   -1.0337];

UpDown = -1; % -1=down, +1=upward

Vel.ens_no = ensemble_number + ensMSB*65536;

% Some will be the same for all ensembles:
Vel.z_adcp =  - UpDown*(dis1/100 + [0:nbins-1]*binlen/100); % meters
Vel.p_adcp = Vel.z_adcp/100; % MPa
Vel.pulselen = pulselen;
Vel.LastByte = LastByte;

% Others will vary:
Vel.yday = datenum(year,month,day,hour,minute,second+hundreths/100) ...
    - datenum(year(1),1,1,0,0,0);

iSV=1:length(Vel.yday);
%iSV=find(Vel.yday>101.7450 & Vel.yday<101.7510); % for PS apr 22,2002 dn/up test
Vel.yday = Vel.yday(iSV);
Vel.ens_no = Vel.ens_no(iSV);
Vel.degC = degC(iSV)/100;
Vel.depth_xducer = xducerdepth(iSV)/10;
Vel.soundvel = soundspeedRDI(iSV);

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
    if isempty(btrange1) % If DN set for NOT bottom-tracking, put in NaNs (3-18-2003)
        ix=[1:length(iSV)];
        Vel.btrange = NaN * [ix; ix; ix; ix];
        Vel.bottomBT = NaN * ix;
        Vel.btvel_bm = NaN * [ix; ix; ix; ix];
        Vel.VbtE = NaN * ix;
        Vel.VbtN = NaN * ix;
        Vel.VbtZ = NaN * ix;
        ANGdn = 90-Vel.SWIMS_headT;
    else % bottom-tracking ON
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
    end
    
end % UpDown<0 for BT

v1=v1(:,iSV);v2=v2(:,iSV);v3=v3(:,iSV);v4=v4(:,iSV);

%% Measured beam components
bad = find(v1==BADVEL); if bad, v1(bad) = NaN; end
bad = find(v2==BADVEL); if bad, v2(bad) = NaN; end
bad = find(v3==BADVEL); if bad, v3(bad) = NaN; end
bad = find(v4==BADVEL); if bad, v4(bad) = NaN; end
Vel.v1_bm = v1/1000; Vel.v2_bm = v2/1000; % m/s
Vel.v3_bm = v3/1000; Vel.v4_bm = v4/1000; % m/s

if savetype>0
    % Save more ...
    Vel.ec1_bm = e1; Vel.ec2_bm = e2; Vel.ec3_bm = e3; Vel.ec4_bm = e4;
    Vel.cor1_bm = cor1; Vel.cor2_bm = cor2; Vel.cor3_bm = cor3; Vel.cor4_bm = cor4;
    if savetype>=2 % implement this later
        vEdn=NaN*Vel.v1_bm; vNdn=vEdn; vZdn=vEdn; vQdn=vEdn;
        for i=1:size(vEdn,1)
            v = TXF*[Vel.v1_bm(i,:); Vel.v2_bm(i,:); Vel.v3_bm(i,:); Vel.v4_bm(i,:)];
            vEdn(i,:) = cos((ANGdn-45)*pi/180).*(v(1,:)) + ...
                sin((ANGdn-45)*pi/180).*(v(2,:));
            vNdn(i,:) = -sin((ANGdn-45)*pi/180).*(v(1,:)) + ...
                cos((ANGdn-45)*pi/180).*(v(2,:));
            vZdn(i,:)=v(3,:); vQdn(i,:)=v(4,:);
        end
        Vel.u_rel = vEdn; Vel.v_rel = vNdn;
        
        %up: flip x(1->2),z directions
        if 0
        vEup=NaN*Vup.v1_bm; vNup=vEup; vZup=vEup; vQup=vEup;
        for i=1:size(vEup,1)
            v = TXF*[Vup.v1_bm(i,:); Vup.v2_bm(i,:); Vup.v3_bm(i,:); Vup.v4_bm(i,:)];
            vEup(i,:) = cos((ANGup-135)*pi/180).*(-v(1,:)) + ...
                sin((ANGup-135)*pi/180).*(v(2,:));
            vNup(i,:) = -sin((ANGup-135)*pi/180).*(-v(1,:)) + ...
                cos((ANGup-135)*pi/180).*(v(2,:));
            vZup(i,:)=-v(3,:); vQup(i,:)=v(4,:);
        end
        end
    end % of "if 0 (savetype>=2)"
end

return



