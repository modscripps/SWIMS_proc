function [CTDsig] = get_SWIMS_SigData(beg_time, end_time, raw_index, raw_path, flag_ok)
% usage: 
%  [CTDsig] = get_SWIMS_SigData(beg_time, end_time, index_file, raw_path);
%  Returns 24-Hz data for a specified time (yearday) range into the 
%   CTDsig structure.  Data are converted to Volts or Hz, but no shifting
%   of scans or calibrations are applied.
%   Structure fields are: Pr,T1Hz,T2Hz,C1Hz,C2Hz,adVolts,status,modCT;
%   timestamps: yday_adj,SWIMStime,SWIMSfasttime;
%   scalars: year, SwimNo.  adVolts is 4 rows for SWIMS1, 8 rows for SWIMS2.
% Inputs are yearday range, Raw-Swims data index file, Path to Raw data;
%
% EG savepath='c:/swims/ps02'; data_dir = fullfile(savepath,'data_mat');
%    IndFld = fullfile(savepath, 'indexes');
%    SIG = get_SWIMS_SigData(118.7, 118.9, ...
%        fullfile(IndFld, 'CTD_ps02_matfiles.mat'), ...
%        fullfile(data_dir, 'CTD') );
% DPW - apr/2002

CTDsig = [];

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

xtra = 15 / 86400; % extra data at ends for sensor lags, filtering (N/A)
SWraw = get_SWIMS_RawData(beg_time-xtra, end_time+xtra, ...
    raw_index, raw_path, flag_ok);

if isempty(SWraw)
    warning('Unable to retrieve raw data')
    return
end

year = SWraw.year;
SwimNo = SWraw.SWIMS_num;
% Use yearday at middle of data interval to determine sensors, calibrations
yday = (SWraw.yday_adj(1)+SWraw.yday_adj(end)) / 2;

%Get Serial numbers.  
[snt1,sparet1]=GetSWIMSConfig('t1',year,yday);
[snt2,sparet2]=GetSWIMSConfig('t2',year,yday);
[snc1,sparec1]=GetSWIMSConfig('c1',year,yday);
[snc2,sparec2]=GetSWIMSConfig('c2',year,yday);
[snpr,spare]=GetSWIMSConfig('pr',year,yday);
[snalt,spare]=GetSWIMSConfig('Altim',year,yday);
% As of 11-mar-04, Dox avail on SWIMS 1, but in Altim channel
[sndox,spare]=GetSWIMSConfig('dox',year,yday); 

if SwimNo==2
    [snflu,spare]=GetSWIMSConfig('flu',year,yday);
    [snobs,spare]=GetSWIMSConfig('obs',year,yday);
    [snpH,spare]=GetSWIMSConfig('pH',year,yday);
    [snH2S,spare]=GetSWIMSConfig('H2S',year,yday);
end

% Save timestamps (in requested range), other info
inRng = find(SWraw.yday_adj>=beg_time & SWraw.yday_adj<=end_time);
CTDsig.yday_adj = SWraw.yday_adj(inRng); % uses elapsed seconds + origin
CTDsig.SWIMStime = SWraw.SWIMStime(inRng); % based on integer seconds
CTDsig.SWIMSfasttime = SWraw.SWIMSfasttime(inRng); % elapsed seconds, float
CTDsig.modCT = SWraw.modCT(inRng);

CTDsig.year = year;
CTDsig.SwimNo = SwimNo;

%Get pressure

CTDsig.Pr = SWraw.Pr; % Already done along with raw counts
CTDsig.status = SWraw.status;

% Check flags for switched C,T cabling (e.g., to activate pumps in fresh water)
if strcmpi(sparet1, 'c1')
    tmp = SWraw.tfreq; SWraw.tfreq = SWraw.cfreq; SWraw.cfreq = tmp;
    %t1_lag=2
end
if strcmpi(sparet2, 'c2')
    tmp = SWraw.tfreq2; SWraw.tfreq2 = SWraw.cfreq2; SWraw.cfreq2 = tmp;
    %t2_lag=2
end
clear tmp

%Get temperature 
CTDsig.T1Hz = SWraw.tfreq/256; 
CTDsig.T2Hz = SWraw.tfreq2/256;
%Get cond 
CTDsig.C1Hz = SWraw.cfreq/256;
CTDsig.C2Hz = SWraw.cfreq2/256;

% A-D data
CTDsig.adVolts = AtoD_SWIMS(SWraw.addata, SwimNo);

% Subsample data to fall in requested yearday range 
vars = {'Pr','T1Hz','T2Hz','C1Hz','C2Hz','adVolts','status'};
for i=1:length(vars)
   % make sure vars{i} exists for this data.
   if ( isfield(CTDsig, vars{i}) )
       eval(['CTDsig.' vars{i} ' = CTDsig.' vars{i} '(:,inRng);'])
   end
end
