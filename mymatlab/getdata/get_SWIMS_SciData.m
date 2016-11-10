function [CTDsci] = get_SWIMS_SciData(beg_time, end_time, raw_index, raw_path, flag_ok, params)
% usage: 
%  [CTDsci] = get_SWIMS_SciData(beg_time, end_time, index_file, raw_path);
%  Returns 24-Hz data for a specified time (yearday) range into the 
%   CTDsci structure.  Calibrations are applied to all active channels:
%   Structure fields are: Pr, T1,T2, C1,C2, Roll,Pitch, [SWIMS2: Dox,Flu,Obs,Obs2(oct/2011)];
%   derived: S1,S2, Th1,Th2, Sg1,Sg2 (salin, theta, sigma_theta),
%   timestamps: yday_adj,SWIMStime,SWIMSfasttime, modCT,
%   scalars: cond_lag, year, SwimNo, several calibration parameters
% Inputs are yearday range, Raw-Swims data index file, Path to Raw data;
%
% EG savepath='c:/swims/ps02'; data_dir = fullfile(savepath,'data_mat');
%    IndFld = fullfile(savepath, 'indexes');
%    SWraw = get_SWIMS_RawData(118.7, 118.9, ...
%        fullfile(IndFld, 'CTD_ps02_matfiles.mat'), ...
%        fullfile(data_dir, 'CTD') );
% DPW - apr/2002; 
% DPW jun-2015: split off processing into function calc_SWIMS_SciData.m

CTDsci = [];

if nargin<3 || isempty(raw_index)
	error(['Inputs beg_time,end_time,raw_index are required.']);
end
if beg_time<0 || end_time>370 || beg_time>=end_time
   error(['beg_time,end_time are out of range/order.']);
end
if nargin<4
   data_path = []; % already in path (??)
end
if nargin<5 || isempty(flag_ok)
   flag_ok = 0;  % don't allow special overrides
end
if nargin<6;
    params = []; % use defaults for lags, calibs, etc.
end
if (end_time-beg_time)*24>12 && flag_ok~=1
   error('Too much data requested, exceeds 12 hours');
end

% check for existence of index file, if it exists load it
% otherwise return with an error and explanation.
if(exist(raw_index) == 0)
   disp('Error:  Index file does not exist.\n');
   error('Specify full or relative path, and include .mat extension');
end

xtra = 15 / 86400; % extra data at ends for sensor lags, filtering
SWraw = get_SWIMS_RawData(beg_time-xtra, end_time+xtra, ...
    raw_index, raw_path, flag_ok);
if isempty(SWraw)
    warning('Unable to retrieve raw data')
    return
end

% process raw data
CTDsci = calc_SWIMS_SciData(SWraw, params);
if isempty(CTDsci)
    return
end

% Subsample data to fall in requested yearday range 
inRng = find(CTDsci.yday_adj>=beg_time & CTDsci.yday_adj<=end_time);
if isempty(inRng)
    CTDsci = [];
    return
end
vars = {'SWIMStime','SWIMSfasttime','yday_adj','modCT','status', ...
    'Pr','T1','T2','C1','C2','S1','S2','Th1','Th2','Sg1','Sg2', ...
    'Roll','Pitch','Dox','Flu','Obs','pH','Sulf','H2S','Alt','Obs2'};
for i=1:length(vars)
   % make sure vars{i} exists for this data.
   if ( isfield(CTDsci, vars{i}) )
       CTDsci.(vars{i}) = CTDsci.(vars{i})(inRng);
   end
end
