function [CTDraw] = get_SWIMS_RawData(beg_time, end_time, index_file, data_path, flag_ok)
% usage: 
%    [CTDraw] = get_SWIMS_RawData(beg_time, end_time, index_file, data_path);
%    returns data for a specified time range into the 
%    CTDraw struct.  Most are raw counts, and date/times, Pressure
%

% example:
%    stuff = get_adcp_data(98, 99, '/a/esp/swims_data/swims/ps01/griddata/', ...
%			'/home/zach/work/ps01/ADCP_Log_ps01.mat', {'u_wat', 'v_wat', 'bottomBT'})
% 
% input:
%    beg_time   - lower limit of the time for which data is to be retrieved. (yeardays)
%    end_time   - upper limit of the time for which data is to be retrieved. (yeardays)
%    index_file - the file which indexes the data files. created by make_ps01_adcp_index.m
%    data_path  - path to the data files which are indexed.
%    vars       - list of the vars of intrest on the data files. ex: { 'u_wat', 'v_wat', 'bottomBT'}
%                 if {} is passed in, the function returns all of the variables.
%
% output:
%    Veldata    - a struct with fields of the requested variables plus 'z_adcp' and 'yday' 
%                 which are always returned.

% revision history:
%    8.20.01    - original version.              zach frazier
%    9.15.01    - bug fixes and modifications.   eileen hash
%    9.18.01    - more bug fixes and header.     zach frazier
%    1.22.02    - Use revised file index, other stuff - Dave W

CTDraw = [];

if nargin<3 | isempty(index_file)
	error(['Inputs beg_time,end_time,index_file are required.']);
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
if (end_time-beg_time)*24>48 & flag_ok~=1
   error('Too much data requested, exceeds 48 hours');
end

vars = []; % make sure all fields are retrieved

% check for existence of index file, if it exists load it
% otherwise return with an error and explanation.
if(exist(index_file) == 0)
   disp('Error:  Index file does not exist.\n');
   error('Specify full or relative path, and include .mat extension');
else
   load(index_file);
end;
iSet=0;
% find first subset in yearday range (cannot straddle subsets)
for i=1:length(Set_params)
   yb = Set_params(i).start_yday;
   ye = Set_params(i).end_yday;
   if yb<end_time & ye>beg_time
      iSet = i;
      break
   end
end
if ~iSet
   error(['Raw Index does not cover specified time range'])
end

% Find and sort files in range
yb = Index(iSet).yday_beg; ye = Index(iSet).yday_end;
ifil = find(yb<end_time & ye>beg_time);
if isempty(ifil)
   warning(['Raw Index has no files for specified time range'])
   return
end
yb = yb(ifil); ye=ye(ifil);
fn=Index(iSet).filename; fn=fn(ifil);
[x,iord] = sort(yb);
fyb = yb(iord); fye = ye(iord); fnam = fn(iord);

% value for	correcting profile-times to actual UTC (EG, sometimes off by 1 hr)
yday_ref = datenum(Set_params(iSet).year,1,1); % start of year
if isfield(Set_params(iSet), 'yday_seconds_offset')
    yd_off = Set_params(iSet).yday_seconds_offset / 86400; 
else
    yd_off = 0;
end

%% Read in SWIMS raw count files
for ifil = 1:length(fyb)
   SWraw = [];
   load( fullfile(data_path, fnam{ifil}) );
   if(isempty(SWraw))
      error(['Error:loading file: ' data_path, fnam{ifil} ]);
   else
      SWraw2(ifil) = SWraw;
   end;
end;

if (length(vars) == 0)
    vars = {'SWIMStime','SWIMSfasttime','tfreq','cfreq','pfreq','tfreq2', ...
            'cfreq2','addata','tprCT','status','modCT','Pr'};
end


% initialize the variables in the new struct.
for i=1:length(vars)
   % first check that vars exist in the structures.
   if( isfield(SWraw, vars{i}) )
       CTDraw = setfield(CTDraw, vars{i},[]);
   else
      disp('Error: vars contains invalid fields.');
      return
   end;
end;

% Return some values for proper calibrations, calculations
CTDraw.year = Set_params(iSet).year;
CTDraw.SWIMS_num = Set_params(iSet).SWIMS_num;

CTDraw.yday_adj = []; % Adjust yeardays using elapsed seconds

for j=1:length(SWraw2)
    clear SWraw
    SWraw = SWraw2(j);
    index = find(SWraw.SWIMStime+yd_off >= beg_time & ...
        SWraw.SWIMStime+yd_off <= end_time & [1 diff(SWraw.SWIMSfasttime)]>0 );
    if ~isempty(index)
        for i=1:length(vars)
            temp = ['[ CTDraw.' vars{i} ', SWraw.' vars{i} '(:,index)];'];
            temp2 = eval(temp);
            CTDraw = setfield(CTDraw, vars{i}, temp2);
        end
        yadj = SWraw.FastTime_Origin + SWraw.SWIMSfasttime(:,index)/86400 + yd_off;
        CTDraw.yday_adj = [CTDraw.yday_adj, yadj];
    end
end;

