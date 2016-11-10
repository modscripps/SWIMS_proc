function Avg_Vel = SWIMS_Vel_Prof(yday_beg, yday_end, cruise) 
% Avg_Vel = SWIMS_Vel_Prof_WaWaves14(yday_beg, yday_end);
% Return results of SWIMS ADCP (UP and DN) averaged over a SWIMS
% profile between yday_beg and yday_end

% Return results of SWIMS ADCP (UP and DN) averaged over a SWIMS
% profile between yday_beg and yday_end.
year = 2016;
adcp_path = '~/Documents/MATLAB/FAST/adcp.mat'; 
    % Path to adcp file from which average profiles will be taken.
depth_range = [20 300]; 
DoPlot = 0;
% End of things you should need to edit!

yd_b = yday_beg; yd_e = yday_end;

try
    AD_SWIMS_PrSynch
catch
    Avg_Vel = []; return
end
    
warning off MATLAB:interp1:NaNinY
try
    AD_SWIMS_SSadj
catch
    Avg_Vel = []; return
end

ADX_clean

ngd = 0; % count how many bins will be usefully averaged
for i=1:length(ADX.depth)
    if length( find( ~isnan(ADX.u(i,:)+ADX.v(i,:)) ) ) > 6
        ngd = ngd + 1;
    end
end
if ngd < 10
    Avg_Vel = []; % not enough bins to bother
    return
end

% Following is for use when SWIMS BT entirely Out-of-range or Off;
% Prepare reference velocity profiles from vessel OS75 ADCP
% (use mean of 5-min VmDas STA's within a few minutes of SWIMS profile)
VM=[];
%    load('/Users/ecfine/Documents/MATLAB/swims/SWIMS_proc/mymatlab/Vels/os75nb_uv.mat') % 5-min os75
    load(adcp_path)
%     VM.yday=sadcp.yday; VM.z_adcp=sadcp.z;
%     VM.u_wat=sadcp.u; VM.v_wat=sadcp.v;
% clear sadcp
    VM.yday=S.datenum-datenum(year,1,1,0,0,0);  VM.z_adcp=S.z;
    VM.u_wat=S.u; VM.v_wat=S.v;
clear S
% VM=get_adcp_any(yd_b-1.6/1440, yd_e+1.1/1440, ...
%      'C:\swims\MORT\sadcp\mat_index\STA_os150_MORT_matfiles.mat', ...
%      'C:\swims\MORT\sadcp\STA_mat');
REF.Z_ok = [0 0]; % depth range to use for seeding SWIMS ADCP avging
if ~isempty(VM)
    if ~isfield(VM, 'yday')
        disp('error coming')
    end
    bot = 0; ig=[];
    ik = find(VM.yday>=yd_b-10/1440 & VM.yday<=yd_e+5/1440);
    if ~isempty(ik)
        REF.Z_ok = [80 111]; % depth range for ship's ADCP ref
        bot = REF.Z_ok(2);
        REF.depth = VM.z_adcp;
        REF.u = nanmean(VM.u_wat(:,ik),2);
        REF.v = nanmean(VM.v_wat(:,ik),2);
        ig = find(~isnan(REF.u+REF.v) & REF.depth<=bot);
    end
    if ~isempty(ig)
        bot = REF.depth(ig(end));
    end
    REF.Z_ok(2) = bot;
end
if isempty(VM) || diff(REF.Z_ok)<10
    REF.depth = [30 40 50];
    REF.Z_ok = [29.999 50.001];
    REF.u = [0; 0; 0]; REF.v = [0; 0; 0];
end
clear VM

% check if bottom-tracking is suitable to determine SWIMS motion
BTokay = 0;  % assume BT off or out-of-range
if length(find(~isnan(ADX.uBT+ADX.vBT+ADX.bottomBT))) > 8
    BTokay = 1; % hope that's enough to use BTing
end
% Compute velocity profiles from SWIMS ADCP cycles
if ~BTokay
    try
        ADX_Prof_avg_ToRef
    catch 
        Avg_Vel = [];
        return 
    end
    AVG.BTokay = 0;
else
    % Use when BT is in range for reasonable number of pings:
    try
        ADX_Prof_avg
        AVG.BTokay = 1;
    catch
        Avg_Vel = [];
        return
    end
end

DoPlot = 0;
x = 0;  % for a break-point to plot in debug mode
if DoPlot
    figure
    subplot(1,2,1)
    plot(REF.u, REF.depth,'k-'), axis ij,grid on,hold on,title('U')
    plot(AVG.U_abs, AVG.depth, 'g-', AVG.U_rel, AVG.depth, 'b-')
    subplot(1,2,2)
    plot(REF.v, REF.depth,'k-'), axis ij,grid on,hold on,title('V')
    plot(AVG.V_abs, AVG.depth, 'g-', AVG.V_rel, AVG.depth, 'b-')
    pause
end

Avg_Vel = AVG;

return