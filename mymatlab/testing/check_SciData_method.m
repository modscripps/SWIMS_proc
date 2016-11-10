% check_SciData_method.m

cruz = 'mc09';
year = 2009; ydb = 105.8; yde = 106.8; % limits for gridded profiles

CalPars = [];

%%% Cruise Defn, folders %%%
savepath = []; % avoid conflict with matlab function
crz=cruz;
set_swims_paths
cruise=cruz;
SWIMSfoldersPC
CrName=cruz;

%IndFld = fullfile(savepath, 'indexes'); % raw,matlab data index folder
IndFld = swimsindex;
%MatFld = fullfile(savepath,'data_mat','CTD'); % local mat directory
MatFld = fullfile(swimsmatdata,'CTD'); % local mat directory
MatIndx = fullfile(IndFld, ['CTD_' CrName '_matfiles.mat']);
GridIndx = fullfile(IndFld, ['SWIMS_' CrName '_gridfiles.mat']);
%GridFld = fullfile(savepath,'griddata'); % local mat directory
GridFld = swimsgridded; % local mat directory

GpsIndx = fullfile(IndFld, ['GPS_' CrName '_matfiles.mat']);
%GpsFld = fullfile(savepath,'data_mat', 'GPS');
GpsFld = fullfile(swimsmatdata,'GPS');

[CTD] = get_swims_data(ydb, yde, ...
    fullfile(swimsindex,['SWIMS_' crz '_gridfiles.mat']), swimsgridded);

return

YLIMS = [];
YLIMS(1,:) = CTD.yday_LDbeg;
YLIMS(2,:) = CTD.yday_LDend;
YLIMS(3,:) = CTD.updown;

yd_b = YLIMS(1,1); yd_e = YLIMS(2,end); % 24-Hz CTD data limits
ZGmin = 0; ZGint = 0.5; ZGmax = 650;
zgrid = [ZGmin:ZGint:ZGmax];
GridOvts = 1;
IgnoreOrder = 1;
IgnoreGPS = 0; 
% get data for gridding
GridStartTime = yd_b;
yd_b = yd_b - (10/86400); % Just before first profile
yd_e = yd_e + (5/86400); % Just after last profile
GPS = get_gps_data(yd_b, yd_e, GpsIndx, GpsFld, 1);
CSci = get_SWIMS_SciData(yd_b, yd_e, MatIndx, MatFld, 1);
    
NGr=[];
GrNum = NaN; % to grid data without affecting index/database

Grid_MakeNewProfs_Calc