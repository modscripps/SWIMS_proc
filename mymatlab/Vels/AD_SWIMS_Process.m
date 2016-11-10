% AD_SWIMS_Process.m

% cruise='home02'; % ADtyp = 'ADDN';
% ADidx = ['C:\swims\' cruise '\indexes\' ADtyp '_' cruise '_matfiles.mat'];
% ADfld = ['C:\swims\' cruise '\data_mat\' ADtyp];
% % yd_b and yd_e should already be defined

ex_yday = 20/86400; % extra data at ends
WC_val = 60-00; %  set beam vel = NaN for correlations below this (default=64)

%% This level of processing makes no corrections for sound speed or clock
%% offsets between the ADCPs and CTD (or GPS).  Time is corrected, however,
%% for jumps in the ADCP clock noted in the matlab-datafile index
%% (i.e., ADDN jumping ahead 1 hour in home02 and bs03).
%% Pitch and roll attitude of the ADCPs is factored in to the velocity
%% computations, and components are saved as:
    %% positive u is flow from NOSE RIGHT to TAIL LEFT (from DN beam 1 to 2),
    %% positive v is flow from TAIL RIGHT to NOSE LEFT (from DN beam 4 to 3),
    %% positive w is flow towards ADDN (SWIMS descending).

ADRaw = get_ADupdn_data(yd_b-ex_yday,yd_e+ex_yday, ADidx, ADfld, 1);

if strcmp(cruise,'wa_nliw_apr2013') && strcmp(ADtyp,'ADUP')
    % special code to copy ADUP pitch and roll from ADDN
    % (ADUP.pitch jumped 4 deg too high from yday=114.0543 to 114.3965)
    AD_SWIMS_Proc_Fix_wanliwapr2013
end

% determine if bottom-track data are present
if isfield(ADRaw, 'btrange')
    BTon = 1;
else
    BTon = 0;
end

% initialize output structure (trim later to exclude extra at ends)
clear ADP
ADP.u = NaN * ADRaw.v1_bm; ADP.v = ADP.u; ADP.w = ADP.u; ADP.werr = ADP.u;
ADP.ec1 = 0*ones(size(ADP.u)); ADP.ec2 = ADP.u; ADP.ec3 = ADP.u; ADP.ec4 = ADP.u;
ADP.yday = ADRaw.yday;
ADP.ens_no = ADRaw.ens_no;
ADP.z_adcp = ADRaw.z_adcp;
if isfield(ADRaw,'depth_xducer')
    ADP.depth_xducer = ADRaw.depth_xducer; % SWIMS 3
else
    ADP.depth_xducer = NaN*ADP.ens_no; % SWIMS 2
end
if size(ADP.z_adcp,1)==1
    ADP.z_adcp = ADP.z_adcp'; % column vector
end

% Screen data for 'gross' problems
ix = find(ADRaw.cor1_bm(:) < WC_val);
ADRaw.v1_bm(ix) = NaN;
ix = find(ADRaw.cor2_bm(:) < WC_val);
ADRaw.v2_bm(ix) = NaN;
ix = find(ADRaw.cor3_bm(:) < WC_val);
ADRaw.v3_bm(ix) = NaN;
ix = find(ADRaw.cor4_bm(:) < WC_val);
ADRaw.v4_bm(ix) = NaN;

% Now set up constants and parameters
theta_o=20; %Angle of the SWIMS ADCP beams from vertical (without pitch or roll)
c_tho = cos(theta_o*pi/180); s_tho = sin(theta_o*pi/180);
SvelA = 1500; % Used in ADCPs' setups
% Make range vector from the depth vector (std, before pitch/roll adjustments)
rBM_o = ADP.z_adcp/c_tho;
% RDI transformation matrix from beams 1-4 to u(1-2), v(4-3), w(avg xz,yz), err vel
Bm2InTx = ...
    [ 1/(2*s_tho) -1/(2*s_tho) 0 0;
      0 0 -1/(2*s_tho) 1/(2*s_tho);
      1/(4*c_tho) 1/(4*c_tho) 1/(4*c_tho) 1/(4*c_tho);
      1/(sqrt(2)*2*s_tho) 1/(sqrt(2)*2*s_tho) -1/(sqrt(2)*2*s_tho) -1/(sqrt(2)*2*s_tho)];

% Turn the recorded pitch into real pitch... insignif. for small rolls.
ADRaw.pitch=180/pi*atan(tan(pi/180*ADRaw.pitch).*cos(pi/180*ADRaw.roll));
% % Slighty filter pitch and roll:
% ADP.pitch = [ ADRaw.pitch(1), ...
%         (0.2*ADRaw.pitch(1:end-2)) + (0.6*ADRaw.pitch(2:end-1)) + (0.2*ADRaw.pitch(3:end)), ...
%         ADRaw.pitch(end) ];
% ADP.roll = [ ADRaw.roll(1), ...
%         (0.2*ADRaw.roll(1:end-2)) + (0.6*ADRaw.roll(2:end-1)) + (0.2*ADRaw.roll(3:end)), ...
%         ADRaw.roll(end) ];

ADP.pitch = ADRaw.pitch; % ensembles sometimes missed, so ...
ADP.roll = ADRaw.roll;  %   just COPY pitch and roll

% Try to eliminate spurious BT returns in each beam
if BTon
    for i=1:4
        mb=medfilt1(ADRaw.btrange(i,:), 5); ix=find(isnan(mb));
        ADRaw.btrange(i,ix) = NaN;
    end
    ADRaw.btrange = ADRaw.btrange/100; % convert to meters
    % initialize output values
    ADP.uBT = NaN * ADRaw.yday;
    ADP.vBT = ADP.uBT; ADP.wBT = ADP.uBT; ADP.werrBT = ADP.uBT;
    ADP.rangeBT = ADP.uBT; ADP.rangeBTmin = ADP.uBT;
end

%% Define values depending on upward or downward looking:
% headings (beam 3) to Magnetic along SWIMS axis (fwd),
clear pD rD
switch ADtyp
    case 'ADDN'
        ADP.SWIMS_headM = mod(ADRaw.heading + (-1)*45, 360);
        sg1 = 1; sg3 = 1;
    case 'ADUP'
        ADP.SWIMS_headM = mod(ADRaw.heading + (1)*45, 360);
        sg1 = 1; sg3 = -1;
end

warning off MATLAB:interp1:NaNinY
for ic=1:length(ADP.yday)
    % The first step is to compute depth vectors for each beam.  These differ from each other 
    % since the beams are tilted.  DO NOT adjust for sound speed (this will happen later).
    % Roll affects beams 1 and 2,  pitch affects beams 3 and 4.
    rBM = rBM_o; % adjust depths for pitched/rolled beam angles
    z1_rel = cos( pi/180*(theta_o + sg1*ADP.roll(ic)) )*rBM;
    z2_rel = cos( pi/180*(theta_o + -sg1*ADP.roll(ic)) )*rBM; 
    z3_rel = cos( pi/180*(theta_o + sg3*ADP.pitch(ic)) )*rBM;
    z4_rel = cos( pi/180*(theta_o + -sg3*ADP.pitch(ic)) )*rBM;
    
    if BTon % correct the bottom ranges for roll and pitch (but not soundspeed yet)
        btrangeCORR(1) = ADRaw.btrange(1,ic) * cos(pi/180*(theta_o+sg1*ADP.roll(ic)))/c_tho;
        btrangeCORR(2) = ADRaw.btrange(2,ic) * cos(pi/180*(theta_o-sg1*ADP.roll(ic)))/c_tho;
        btrangeCORR(3) = ADRaw.btrange(3,ic) * cos(pi/180*(theta_o+sg3*ADP.pitch(ic)))/c_tho;
        btrangeCORR(4) = ADRaw.btrange(4,ic) * cos(pi/180*(theta_o-sg3*ADP.pitch(ic)))/c_tho;
        %Average all four together to get our best estimate of bottom depth below SWIMS
        ig = find(~isnan(btrangeCORR));
        if length(ig)>1
            ADP.rangeBT(ic) = mean(btrangeCORR(ig));
            ADP.rangeBTmin(ic) = min(btrangeCORR(ig));
        end
    end
    
    % map beam velocities onto standard depth-grid (extend first bin up, if needed)
    u1 = interp1( [0; z1_rel], [ADRaw.v1_bm(1,ic); ADRaw.v1_bm(:,ic)], ADP.z_adcp);
    u2 = interp1( [0; z2_rel], [ADRaw.v2_bm(1,ic); ADRaw.v2_bm(:,ic)], ADP.z_adcp);
    u3 = interp1( [0; z3_rel], [ADRaw.v3_bm(1,ic); ADRaw.v3_bm(:,ic)], ADP.z_adcp);
    u4 = interp1( [0; z4_rel], [ADRaw.v4_bm(1,ic); ADRaw.v4_bm(:,ic)], ADP.z_adcp);
    % map echo intensities, save for later screening (extend first bin up, if needed)
    ADP.ec1(:,ic) = round( interp1( [0; z1_rel], [ADRaw.ec1_bm(1,ic); ADRaw.ec1_bm(:,ic)], ...
        ADP.z_adcp, 'linear', 0) );
    ADP.ec2(:,ic) = round( interp1( [0; z2_rel], [ADRaw.ec2_bm(1,ic); ADRaw.ec2_bm(:,ic)], ...
        ADP.z_adcp, 'linear', 0) );
    ADP.ec3(:,ic) = round( interp1( [0; z3_rel], [ADRaw.ec3_bm(1,ic); ADRaw.ec3_bm(:,ic)], ...
        ADP.z_adcp, 'linear', 0) );
    ADP.ec4(:,ic) = round( interp1( [0; z4_rel], [ADRaw.ec4_bm(1,ic); ADRaw.ec4_bm(:,ic)], ...
        ADP.z_adcp, 'linear', 0) );
    
    % Now, transform beam velocities to RDI instrument frame components.
    % Manipulate ADUP beams, pitch, roll to match ADDN frame definition:
    %   UP beams 1,2,3,4 are same directions as DN beams -4,-3,-2,-1 and
    %   UP pitch = DN roll, UP roll = -DN pitch (for SWIMS setup).
    % First, into horizontal, vertical components in tilted ADDN-frame,
    % then into leveled components.
    switch ADtyp
        case 'ADDN'
            VelsIN = ( Bm2InTx * [u1 u2 u3 u4]' )';
            pD = ADP.pitch(ic); rD = ADP.roll(ic);
        case 'ADUP'
            VelsIN = ( Bm2InTx * [-u4 -u3 -u2 -u1]' )'; % to match ADDN
            pD = -ADP.roll(ic); rD = ADP.pitch(ic);
    end
    % Rotate by (-roll) about 4-3 DN-axis, then by (-pitch) about leveled 1-2 DN-axis:
    In2Lev = ...
        [ 1 0 0;
          0 cos(pi/180*pD) -sin(pi/180*pD);
          0 sin(pi/180*pD) cos(pi/180*pD) ] * ...
        [ cos(pi/180*rD) 0 sin(pi/180*rD);
          0 1 0;
          -sin(pi/180*rD) 0 cos(pi/180*rD) ];
    VelsLev = ( In2Lev * VelsIN(:,1:3)' )';
    % Gather components into output structure
    ADP.u(:,ic) = VelsLev(:,1);
    ADP.v(:,ic) = VelsLev(:,2);
    ADP.w(:,ic) = VelsLev(:,3);
    ADP.werr(:,ic) = VelsIN(:,4);

    % Also for bottom track velocity:
    if BTon
        VelsBTIN = ( Bm2InTx * ADRaw.btvel_bm(1:4,ic) )';
        VelsBTLev = ( In2Lev * VelsBTIN(:,1:3)' )';
        ADP.uBT(:,ic) = VelsBTLev(:,1);
        ADP.vBT(:,ic) = VelsBTLev(:,2);
        ADP.wBT(:,ic) = VelsBTLev(:,3);
        ADP.werrBT(:,ic) = VelsBTIN(:,4);
    end
    
end

clear ADRaw
%% Trim to original yearday range
vn = fieldnames(ADP);
id = find( ADP.yday>=yd_b & ADP.yday<=yd_e );
for iv=1:length(vn)
    if ~strcmp(vn{iv}, 'z_adcp')
        ADP.(vn{iv}) = ADP.(vn{iv})(:,id);
    end
end