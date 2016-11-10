% set_swims_paths.m

% -------------------------------------------------------------------------
% Edit for your personal setup and cruise
% -------------------------------------------------------------------------

crz = 'SproulTest';
year = 2016;
%directory containing raw data
 %localpath   = fullfile('/Users','ecfine','Documents','scripps',...
 %    'arctic2015','swims','raw');
 localpath   = fullfile('/Users','ecfine', 'Documents','FAST','raw_data');
   
matlabdisk  = fullfile('/Users','ecfine','Documents','MATLAB');
%directory where processed data is saved 
savepath=fullfile(matlabdisk,'swims', crz);

% Convenient to have for loading data
index_file  = sprintf('%s/indexes/SWIMS_%s_gridfiles.mat',savepath,crz);
data_path   = fullfile(savepath,'griddata');
vars        = {'p','th1','th2','s1','s2','sgth1','sgth2','eps2','yday',...
    'lat','lon','t1','t2','c1','c2','flu','dox'};


% -------------------------------------------------------------------------
% Nothing below this should be altered!
% -------------------------------------------------------------------------

mymatlab    = fullfile(matlabdisk,'swims','SWIMS_proc','mymatlab'); 

swimsindex=fullfile(matlabdisk, 'swims', crz, 'indexes');
swimsthalweg=fullfile(matlabdisk, 'swims', crz, 'thalwegs');
swimsgridded=fullfile(matlabdisk, 'swims', crz, 'griddata');
swimsmatdata=fullfile(matlabdisk, 'swims', crz, 'data_mat');
cruise = crz;

global SWIMS_cal
SWIMS_cal=fullfile(matlabdisk,'swims','SWIMS_proc','SWIMS_config','SWIMS_cal');
global SWIMS_config
SWIMS_config=fullfile(matlabdisk,'swims','SWIMS_proc','SWIMS_config','SWIMS_config');
remotepath = localpath;

path(path,matlabdisk)
path(path,mymatlab)
path(fullfile(matlabdisk,'toolboxes'),path)
path(fullfile(matlabdisk,'toolboxes','seawater'),path)
path(fullfile(mymatlab,'Calibration'),path)
path(fullfile(mymatlab,'FieldApp'),path)
path(fullfile(mymatlab,'ReadRoutines'),path)
path(fullfile(mymatlab,'Utilities'),path)
path(fullfile(mymatlab,'getdata'),path)
path(fullfile(mymatlab,'testing'),path)
path(fullfile(mymatlab,'read_raw'),path)
path(fullfile(mymatlab,'overturns'),path)
path(fullfile(mymatlab,'MHAUtilities'),path)
path(fullfile(mymatlab,'FieldApp','SWIMSFieldApp2'),path)
path(fullfile(mymatlab,'FieldApp','SWIMSFieldApp2','MHAUtils'),path)
path(fullfile(mymatlab,'FieldApp','SWIMSFieldApp2','MHA_PS'),path)
path(fullfile(mymatlab,'FieldApp','SWIMSFieldApp2','PublishingTools'),path)
path(fullfile(matlabdisk,'swims','plot'),path)
