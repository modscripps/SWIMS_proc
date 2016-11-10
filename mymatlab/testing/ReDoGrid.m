% ReDoGrid.m - recompute some gridded data, no effect on grid database

clear

crz='stf07'; year = 2007; %'ML04'; %'hc03'; % crz='ps03';
set_swims_paths
cruise='stf07'; %'ML04'; %'hc03'; % cruise='ps03';
SWIMSfoldersPC
CrName='stf07'; %'ML04'; %'hc03'; % CrName = 'ps03'; %'home02'; % 'ps02';

IndFld = fullfile(savepath, 'indexes'); % raw,matlab data index folder
MatFld = fullfile(savepath,'data_mat','CTD'); % local mat directory
MatIndx = fullfile(IndFld, ['CTD_' CrName '_matfiles.mat']);
GridIndx = fullfile(IndFld, ['SWIMS_' CrName '_gridfiles.mat']);
GridFld = fullfile(savepath,'griddata'); % local mat directory

% check portion with S1 all NaNs - 07/07/07
yd_b = 187.128; yd_e = 187.2014;
load(fullfile(GridFld,'SWIMSgrid-2007-187-0305.mat'))

iP = find(SWIMSgrid.yday>=yd_b & SWIMSgrid.yday<yd_e);
nSamps = SWIMSgrid.ends(iP) - SWIMSgrid.starts(iP) + 1;
ydbs = SWIMSgrid.yday(iP);
uds = SWIMSgrid.updown(iP);
ydes = NaN*ydbs;
zgrid = SWIMSgrid.z;
clear SWIMSgrid

yd_b = yd_b - (60/86400); % Just before first profile
yd_e = yd_e + (60/86400); % Just after last profile

CSci = get_SWIMS_SciData(yd_b, yd_e, MatIndx, MatFld, 1);

for i=1:length(iP)
    iy = find(CSci.yday_adj>=ydbs(i)); iy = iy(1);
    ydes(i) = CSci.yday_adj( iy+nSamps(i)-1 );
end

YLIMS = [ydbs; ydes; uds];
GrNum = NaN; % so that there is no effect on the Grid database

Grid_MakeNewProfs_Calc
% NGr will be the freshly gridded stuff
