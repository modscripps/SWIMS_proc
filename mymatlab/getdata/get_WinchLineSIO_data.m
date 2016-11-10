function [Ldata] = get_WinchLineSIO_data(beg_time, end_time, index_file, data_path, flag_ok)
% usage: 
%    [Ldata] = get_WinchLineSIO_data(beg_time, end_time, index_file, data_path);
%    returns line-out, line-rate data for a specified time range into the 
%    Ldata struct.
%
% input:
%    beg_time   - lower limit of the time for which data is to be retrieved. (yeardays=sattime)
%    end_time   - upper limit of the time for which data is to be retrieved. (yeardays=sattime)
%    index_file - the file which indexes the data files.
%    data_path  - path to the data files which are indexed.
%

Ldata = [];

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
if (end_time-beg_time)*24>96 && flag_ok~=1
   error('Too much data requested, exceeds 48 hours');
end

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
   disp(['get_WinchLineSIO_data: Index does not cover specified time range'])
   return
end

%%%% from adcp code: %%%
% Two types of time corrections: yd_off is gross adjustment for mis-set PC clock;
%  yd_jump synchronizing non-UTC timestamps on serial data (LD here) to
%  UTC timestamps in the CTD data (su.exe assigns these differently, 6/07)
% This is NOT intended to account for differing clock drifts of PC vs ADCPs vs GPS.
% NOTE: when this happens, the clock should not be reset until SWIMS will be out of
%      the water for (e.g.) yd_jump or more (to avoid ensemble time conflicts)
yd_jump_aft = [];
yd_jump_by = [];
if isfield(Set_params(iSet), 'yday_jump_after')
    yd_jump_aft = [Set_params(iSet).yday_jump_after,  inf]; % last one to adjust to end
    yd_jump_by = Set_params(iSet).yday_jump_adjust;
end
yd_off = [];
if isfield(Set_params(iSet), 'ping_seconds_offset')
    yd_off = Set_params(iSet).ping_seconds_offset / 86400; 
elseif isfield(Set_params(iSet), 'yday_seconds_offset')
    yd_off = Set_params(iSet).yday_seconds_offset / 86400; 
end
if isempty(yd_off)
    yd_off = 0;
end

% First adjust the time jumps (based on recorded ADCP time), then apply general offset

% Find and sort files in range (adjust times in index first)
yb = Index(iSet).yday_beg; ye = Index(iSet).yday_end;
for i=1:length(yd_jump_by)
    ix = find( yb>=yd_jump_aft(i) & yb<yd_jump_aft(i+1) );
    yb(ix) = yb(ix) + yd_jump_by(i);
    ix = find( ye>=yd_jump_aft(i) & ye<yd_jump_aft(i+1) );
    ye(ix) = ye(ix) + yd_jump_by(i);
end
yb = yb + yd_off; ye = ye + yd_off;

ifil = find(yb<end_time & ye>beg_time);
if isempty(ifil)
   disp(['get_WinchLine_data: Index has no files for specified time range'])
   return
end

yb = yb(ifil); ye=ye(ifil);
fn=Index(iSet).filename; fn=fn(ifil);
[x,iord] = sort(yb);
fyb = yb(iord); fye = ye(iord); fnam = fn(iord);
fn = [];

yday_ref = datenum(Set_params(iSet).year,1,1); % start of year

%% Read in LD files, correct for midnight problem 
for ifil = 1:length(fyb)
    LD = [];
    load( fullfile(data_path, fnam{ifil}) );
    if(isempty(LD))
        error(['Error:loading file: ' data_path, fnam{ifil} ]);
    else
        Ld2(ifil) = LD;
    end;
end;

vars = {'PCyday', 'LCIyday', 'line_rate', 'line_out', 'TensionLbs'};

% initialize the variables in the new struct.
for i=1:length(vars)
   % first check that vars exist in the structures.
   if( isfield(LD, vars{i}) )
       Ldata = setfield(Ldata, vars{i},[]);
   else
      error('Error: File is missing expected field.');
   end;
end;
% Gather data from files into output structure

for j=1:length(Ld2)
    clear LD
    LD = Ld2(j);
    % Adjust timestamps WITHIN data
    for i=1:length(yd_jump_by)
        ix = find( LD.PCyday>=yd_jump_aft(i) & ...
            LD.PCyday<yd_jump_aft(i+1) );
        LD.PCyday(ix) = LD.PCyday(ix) + yd_jump_by(i);
    end
    LD.PCyday = LD.PCyday + yd_off;
    % save data (adjusted) in yearday range
    index = find(LD.PCyday >= beg_time & LD.PCyday <= end_time);
    if ~isempty(index)
        for i=1:5
            temp = ['[ Ldata.' vars{i} ', LD.' vars{i} '(:,index)]'];
            temp2 = eval(temp);
            Ldata = setfield(Ldata, vars{i}, temp2);
        end;
    end

end;
