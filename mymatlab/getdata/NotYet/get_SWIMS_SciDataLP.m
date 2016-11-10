function [CTDsci] = get_SWIMS_SciDataLP(beg_time, end_time, raw_index, raw_path, flag_ok)
% usage: 
%  [CTDsci] = get_SWIMS_SciDataLP(beg_time, end_time, index_file, raw_path);
%  Returns 24-Hz data for a specified time (yearday) range into the 
%   CTDsci structure.  Calibrations are applied to all active channels:
%   Structure fields are: Pr, T1,T2, C1,C2, Roll,Pitch, [SWIMS2: Dox,Flu,Obs];
%   derived: S1,S2, Th1,Th2, Sg1,Sg2 (salin, theta, sigma_theta),
%   timestamps: yday_adj,SWIMStime,SWIMSfasttime, modCT,
%   scalars: cond_lag, year, SwimNo.
%  Some quantities are then low-pass filtered (only Pr for now, mar-2003)
% Inputs are yearday range, Raw-Swims data index file, Path to Raw data;
%
% EG savepath='c:/swims/BS03'; data_dir = fullfile(savepath,'data_mat');
%    IndFld = fullfile(savepath, 'indexes');
%    SWraw = get_SWIMS_SciDataLP(81.1, 81.2, ...
%        fullfile(IndFld, 'CTD_bso3_matfiles.mat'), ...
%        fullfile(data_dir, 'CTD') );
% DPW - Mar/2003

CTDsci = [];

if nargin<3 | isempty(raw_index)
	error(['Inputs beg_time,end_time,raw_index are required.']);
end
if beg_time<0 | end_time>370 | beg_time>=end_time
   error(['beg_time,end_time are out of range/order.']);
end
if nargin<4
   data_path = []; % already in path (??)
end
if nargin<5
   flag_ok = 0;  % don't allow special overrides
end
if (end_time-beg_time)*24>12 & flag_ok~=1
   error('Too much data requested, exceeds 12 hours');
end

% check for existence of index file, if it exists load it
% otherwise return with an error and explanation.
if(exist(raw_index) == 0)
   disp('Error:  Index file does not exist.\n');
   error('Specify full or relative path, and include .mat extension');
end

xtra = 15 / 86400; % extra data at ends for sensor lags, filtering
SD = get_SWIMS_SciData(beg_time-xtra, end_time+xtra, ...
    raw_index, raw_path, flag_ok)
if isempty(SD)
    error('Unable to retrieve unfiltered data')
end

SecBW = 2; % length of 4th-order Butterworth filter to lowpass pressure data
SecWCD = 2; % length in seconds for computing center-difn fallrates

HzSamp = 24; % sampling frequency of pres,yday data

PtoZ = 100; % factor to convert pressures to depth in meters
MetProf = 5; % minimum length (m) for profiles
SecGap = 3; % sample gaps (recorded) exceeding this are flagged as discontinuities
SecTSer = 8; % minimum length (s) for time series
WmsProf = 0.05; % minimum profiling speed (m/s)

% Compute (even) number of samples >= filter times
NumBW = ceil(SecBW*HzSamp/2)*2; 
NumWCD = ceil(SecWCD*HzSamp/2)*2;
NhfBW = NumBW/2; NhfWCD = NumWCD/2; % half lengths
% Get filter coefs
if NumBW >= 2
    [bBW,aBW] = MHAButter(1/HzSamp ,SecBW);
end

% Force column vectors, exclude bad pressure or time data
if size(pr,1) == 1
    pr = pr';
end
if size(yday,1) == 1
    yday = yday';

year = SWraw.year;
SwimNo = SWraw.SWIMS_num;
yday = SWraw.yday_adj(1);


% Save timestamps (in requested range), other info
inRng = find(SWraw.yday_adj>=beg_time & SWraw.yday_adj<=end_time);
CTDsci.yday_adj = SWraw.yday_adj(inRng); % uses elapsed seconds + origin
CTDsci.SWIMStime = SWraw.SWIMStime(inRng); % based on integer seconds
CTDsci.SWIMSfasttime = SWraw.SWIMSfasttime(inRng); % elapsed seconds, float
CTDsci.modCT = SWraw.modCT(inRng);

CTDsci.cond_lag = cond_lag;
CTDsci.dox_lag = dox_lag;
CTDsci.flu_lag = flu_lag;
CTDsci.pH_lag = pH_lag;
CTDsci.h2s_lag = h2s_lag;
CTDsci.year = year;
CTDsci.SwimNo = SwimNo;

%Get pressure

CTDsci.Pr = SWraw.Pr; % Already done along with raw counts

% Subsample data to fall in requested yearday range 
vars = {'Pr','T1','T2','C1','C2','S1','S2','Th1','Th2','Sg1','Sg2', ...
        'Roll','Pitch','Dox','Flu','Obs','pH','Sulf','H2S'};
for i=1:length(vars)
   % make sure vars{i} exists for this data.
   if ( isfield(CTDsci, vars{i}) )
       eval(['CTDsci.' vars{i} ' = CTDsci.' vars{i} '(inRng);'])
   end
end
