function [GPSdata] = get_gps_data(beg_time, end_time, index_file, data_path, flag_ok)
% usage: 
%    [GPSdata] = get_gps_data(beg_time, end_time, index_file, data_path);
%    returns data from the ps01 cruise for a specified time range into the 
%    GPSdata struct.
%
% example:
%    stuff = get_gps_data(98, 99, '/a/esp/swims_data/swims/ps01/griddata/', ...
%			'/home/zach/work/ps01/GPS_Log_ps01.mat')
% 
% input:
%    beg_time   - lower limit of the time for which data is to be retrieved. (yeardays=sattime)
%    end_time   - upper limit of the time for which data is to be retrieved. (yeardays=sattime)
%    index_file - the file which indexes the data files. created by make_ps01_gps_index.m
%    data_path  - path to the data files which are indexed.
%
% output:
%    GPSdata    - a struct with fields of 'lat', 'lon', and 'sattime'

% revision history:
%    8.17.01    - original version.              zach frazier
%    1.22.02    - Use revised file index, other stuff - Dave W

GPSdata = [];

if nargin<3 | isempty(index_file)
	disp(['get_gps_data: Inputs beg_time,end_time,index_file are required.'])
    return
end
if beg_time<0 | end_time>370 | beg_time>=end_time
   disp(['get_gps_data: beg_time,end_time are out of range/order.'])
   return
end
if nargin<4
   data_path = []; % already in path (??)
end
if nargin<5
   flag_ok = 0;  % don't allow special overrides
end
if (end_time-beg_time)*24>96 & flag_ok~=1
   disp('get_gps_data: Too much data requested, exceeds 48 hours');
   return
end

% check for existence of index file, if it exists load it
% otherwise return with an error and explanation.
if(exist(index_file) == 0)
   disp('get_gps_data:  Index file does not exist.\n');
   disp('Specify full or relative path, and include .mat extension');
   return
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
   disp(['get_gps_data:  Index does not cover specified time range'])
   return
end

% Find and sort files in range
yb = Index(iSet).yday_beg; ye = Index(iSet).yday_end;
ifil = find(yb<end_time & ye>beg_time);
if isempty(ifil)
   warning(['Index has no files for specified time range'])
   return
end
yb = yb(ifil); ye=ye(ifil);
fn=Index(iSet).filename; fn=fn(ifil);
[x,iord] = sort(yb);
fyb = yb(iord); fye = ye(iord); fnam = fn(iord);

yday_ref = datenum(Set_params(iSet).year,1,1); % start of year

%% Read in GPS files, correct for midnight problem 
for ifil = 1:length(fyb)
    GPS = [];
    load( fullfile(data_path, fnam{ifil}) );
    if(isempty(GPS))
        error(['Error:loading file: ' fullfile(data_path, fnam{ifil}) ]);
    else
        ix = find(GPS.sattime<GPS.sattime(1));
        GPS.sattime(ix) = GPS.sattime(ix)+1;
        Gps2(ifil) = GPS;
    end;
end;

vars = {'lat','lon','PCyday'};

% initialize the variables in the new struct.
for i=1:length(vars)
   % first check that vars exist in the structures.
   if( isfield(GPS, vars{i}) )
       GPSdata = setfield(GPSdata, vars{i},[]);
   else
      disp('get_gps_data: File is missing expected field.');
      GPSdata = [];
      return
   end;
end;
% Gather data from files into output structure
GPSdata.sattime = [];
for j=1:length(Gps2)
    clear GPS
    GPS = Gps2(j);
    index = find(GPS.sattime >= beg_time & GPS.sattime <= end_time);
    if ~isempty(index)
        for i=1:length(vars)
            temp = ['[ GPSdata.' vars{i} ', GPS.' vars{i} '(:,index)]'];
            temp2 = eval(temp);
            GPSdata = setfield(GPSdata, vars{i}, temp2);
        end;
        GPSdata.sattime = [GPSdata.sattime  GPS.sattime(index)];
    end
end;
