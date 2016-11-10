% Run_SWIMS_VelProfs_byGrid.m

% Initialize structure for saving SWIMS Velocity profiles and related info,
% Compute from pre-processed VelDN,UP files for each profile in a grid file,
% then store shear results with velocity results in corresponding file.
ADProfVars = { 'U_abs','V_abs','U_rel','V_rel',...
        'U_abs_std','V_abs_std','U_rel_std','V_rel_std',...
        'count_abs','count_rel','outlier_abs','outlier_rel',...
        'U_rel2abs','V_rel2abs', 'yd_min','yd_max','BTokay',...
        'dUdz','dVdz','dU_std','dV_std','count_dUV','outlier_dUV'};

set_swims_paths
IX = fullfile(swimsindex,['SWIMS_' crz '_gridfiles.mat']);
DG = swimsgridded;
VG = fullfile(swimsmatdata,'VelSW');
clear SWIMSgrid
MIN_SEC = 12; % minimum duration(sec) of up/downcast to compute vels

load(IX)
lGyd = 0;
lGprYD = -0; % 0 to start, then set = -(last 1-Hz yday) to skip to newer grid profiles
iSet = 1; % SWIMS 2 period
iN0 = 1;
for iN=iN0:length(Index(iSet).filename)
    x = Index(iSet).filename{iN};
    if ~strcmp(x,  'SWIMSgrid-20150828T000258.mat')...
            ...%&&~strcmp(x, 'SWIMSgrid-20150905T035553.mat')...
            ...%&&~strcmp(x,'SWIMSgrid-20150905T075134.mat')...
            &&~strcmp(x,'SWIMSgrid-20150912T000445.mat')...
            &&~strcmp(x,'SWIMSgrid-20150912T053958.mat')...
            &&~strcmp(x,'SWIMSgrid-20150912T155239.mat')...
            &&~strcmp(x,'SWIMSgrid-20150914T012541.mat')...
            &&~strcmp(x,'SWIMSgrid-20150915T093259.mat') %"the grid vectors are not strictly monotonic inc"
%          &&~strcmp(x,'SWIMSgrid-20160317T183829.mat')...
%              &&~strcmp(x,'SWIMSgrid-20160317T211025.mat')...
%               &&~strcmp(x,'SWIMSgrid-20160318T015212.mat')...
%                &&~strcmp(x,'SWIMSgrid-20160318T061247.mat')...
%                 &&~strcmp(x,'SWIMSgrid-20160318T115405.mat')...
%                  &&~strcmp(x,'SWIMSgrid-20160318T155645.mat')...
%                   &&~strcmp(x,'SWIMSgrid-20160318T200025.mat')...
%                    &&~strcmp(x,'SWIMSgrid-20160318T235916.mat')...
%                     &&~strcmp(x,'SWIMSgrid-20160319T120222.mat')...
%                     &&~strcmp(x,'SWIMSgrid-20160319T160214.mat')...

        
        
        load(fullfile(DG, x))
        disp(['Velocity profiles for ' x ' ...'])
        x = ['VelSW' x(10:end)];
        Vfn = fullfile(VG, x); % save velocities here

        ct = length(SWIMSgrid.yday_LDend);
        yd_begs = SWIMSgrid.yday_LDbeg;
        yd_ends = SWIMSgrid.yday_LDend;

        % Prepare structures
        clear SWIMS_Vels
        for i=1:length(ADProfVars)
            SWIMS_Vels.(ADProfVars{i}) = [];
        end
        SWIMS_Vels.depth = [];
        SWIMS_Vels.depSH = [];
        SWIMS_Vels.yday_beg = [];
        SWIMS_Vels.yday_end = [];
        SWIMS_Vels.yday = []; % same as yday_beg, for get_ADCP_any()
        clear SWIMSgrid

        for ip = 1:ct
            if (yd_ends(ip)-yd_begs(ip))*86400 < MIN_SEC || isnan(yd_ends(ip)-yd_begs(ip))
                disp([ip yd_begs(ip) yd_ends(ip)])
                continue % not long enough
            end
            yd_b = yd_begs(ip) - 5/86400; % some extra (5 sec) at ends
            yd_e = yd_ends(ip) + 5/86400;
            clear AVel
            AVel = SWIMS_Vel_Prof(yd_b, yd_e, crz);
            if isempty(AVel)
                disp([-ip yd_begs(ip) yd_ends(ip)])
                continue; % not enough data to average
            end
            pause(0.2)
            % Save ancillary info, regarding ensembles going into averages
            if isempty(SWIMS_Vels.depth)
                SWIMS_Vels.depth = AVel.depth;
                SWIMS_Vels.depSH = AVel.depSH;
                SWIMS_Vels.Info = AVel.ENSinfo;
            else
                n = length(SWIMS_Vels.Info) + 1;
                SWIMS_Vels.Info(n) = AVel.ENSinfo;
            end
            %% accumulate average profiles and related stats
            for i=1:length(ADProfVars)
                SWIMS_Vels.(ADProfVars{i}) = ...
                    [SWIMS_Vels.(ADProfVars{i}) AVel.(ADProfVars{i})];
            end
            % yearday span, also (without extra at ends)
            SWIMS_Vels.yday_beg(end+1) = yd_begs(ip);
            SWIMS_Vels.yday_end(end+1) = yd_ends(ip);
            SWIMS_Vels.yday(end+1) = yd_begs(ip);

            clear AVel
        end
        save(Vfn, 'SWIMS_Vels')
        clear SWIMS_Vels
    end
end

% THEN, make INDEX with script:
%   fullfile(swimsmatdata,'VelSW','Update_VelSW_index.m')

