function [ADdata] = get_ADupdn_data(beg_time, end_time, index_file, data_path, flag_ok)
%    [ADdata] = get_ADupdn_data(beg_time, end_time, index_file, data_path);
%    returns data from SWIMS-mounted ADCPs for a specified time range into the 
%    ADdata struct.  Profile data are NOT included.
% example:
%    AD = get_ADupdn_data(124.8, 124.9, ...
%       'C:\swims\ps02\indexes\ADDN_ps02_matfiles.mat', ...
%		'C:\swims\ps02\data_mat\ADDN\');
% input:
%    beg_time   - lower yearday limit for which data are to be retrieved.
%    end_time   - upper yearday limit for which data are to be retrieved.
%    index_file - the file which indexes the matlab data files.
%    data_path  - path to the matlab data files which are indexed.
% output:
%    ADdata    - a struct with fields of 'ens_no', 'yday', and pitch,roll,
%   and heading info, both in ADCP frame and in SWIM's. ADDN also has
%   bottom track ranges and velocities (not corrected for pitch/roll).

% revision history:
%    20-sep-2003    - original operational version.              Dave W

ADdata = [];

if nargin<3 | isempty(index_file)
	disp(['Inputs beg_time,end_time,index_file are required.']);
    return
end
if beg_time<0 | end_time>370 | beg_time>=end_time
   disp(['beg_time,end_time are out of range/order.']);
   return
end
if nargin<4
   data_path = []; % already in path (??)
end
if nargin<5
   flag_ok = 0;  % don't allow special overrides
end
if (end_time-beg_time)*24>48 & flag_ok~=1
   disp('Too much data requested, exceeds 48 hours');
   return
end

% check for existence of index file, if it exists load it
% otherwise return with an error and explanation.
if(exist(index_file) == 0)
   disp('Error:  Index file does not exist.');
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
   disp(['Index does not cover specified time range'])
   return
end

% Two types of time corrections: yd_off is gross adjustment for mis-set ADCP clock;
%  yd_jump accounts for the ADDN clock suddenly jumping ahead 1 hour from one ping
%      to the next (and staying that way until manually reset).
% This is NOT intended to account for differing clock drifts of PC vs ADCPs vs GPS.
% NOTE: when this happens, the clock should not be reset until SWIMS will be out of
%      the water for (e.g.) 1 hour or more (to avoid ensemble time conflicts)
yd_jump_aft = [];
yd_jump_by = [];
if isfield(Set_params(iSet), 'yday_jump_after')
    yd_jump_aft = [Set_params(iSet).yday_jump_after,  inf]; % last one to adjust to end
    yd_jump_by = Set_params(iSet).yday_jump_adjust;
end
yd_off = Set_params(iSet).yday_seconds_offset / 86400;
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
   disp(['Index has no files for specified time range'])
   return
end

yb = yb(ifil); ye=ye(ifil);
fn=Index(iSet).filename; fn=fn(ifil);
[x,iord] = sort(yb);
fyb = yb(iord); fye = ye(iord); fnam = fn(iord);

yday_ref = datenum(Set_params(iSet).year,1,1); % start of year

%% Read in ADCP files 
clear AD, AD = [];
for ifil = 1:length(fyb)
    SN = []; XX = [];
    XX = load( fullfile(data_path, fnam{ifil}) );
    if ~isempty(XX), SN=fieldnames(XX); end
    if ~isempty(SN) & ~( strcmp(SN{1},'ADUP') | strcmp(SN{1},'ADDN') | ...
            strcmp(SN{1},'VelUP') | strcmp(SN{1},'VelDN') )
        disp(['Error: loading file: ' fullfile(data_path, fnam{ifil}) ]);
        return
    else
        SN=SN{1};
        eval(['fn = fieldnames(XX.' SN ');']);
        if ifil==1
            eval(['AD = XX.' SN ';'])
            fn0 = fn; % keep intersection of fields
        else
            fn0c = fn0;
            for i=1:length(fn0)
                ok = eval(['isfield(XX.' SN ', ''' fn0{i} ''')']);
                if ok
                    if ~strcmp(fn0{i}, 'z_adcp')
                        eval(['AD.' fn0{i} ' = [AD.' fn0{i} ', XX.' SN '.' fn0{i} '];'])
                    end
                else
                    eval(['AD.' fn0{i} ' = [];'])
                    fn0c(i) = []; % remove from further accumulating
                    warning(['Excluding field = ' fn0{i} ' not in all ' SN ' files.'])
                end
            end
            fn0 = fn0c; % update if any removed
        end % 
    end % of concatinating fields from current file
    clear XX
end % of all files in yearday range

if ~isstruct(AD)
    return
end

% Adjust ensemble times WITHIN data
for i=1:length(yd_jump_by)
    ix = find( AD.yday>=yd_jump_aft(i) & AD.yday<yd_jump_aft(i+1) );
    AD.yday(ix) = AD.yday(ix) + yd_jump_by(i);
end
AD.yday = AD.yday + yd_off;

vars = fieldnames(AD);
vars = vars(find(~strcmp(vars, 'z_adcp')));
% ADdata.z_adcp = AD.z_adcp;

index = find(AD.yday >= beg_time & AD.yday <= end_time);
if ~isempty(index)
    ADdata.z_adcp = AD.z_adcp;
    for i=1:length(vars)
        if ~isempty(vars{i})
            ok = eval(['size(AD.' vars{i} ', 2) == length(AD.yday)']);
            if ok  % right size to save and return
                temp = ['ADdata.' vars{i} ' = AD.' vars{i} '(:,index);'];
                eval(temp);
            end
        end
    end
end


