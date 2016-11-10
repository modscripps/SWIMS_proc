function [Veldata] = get_adcp_data(beg_time, end_time, index_file, data_path, vars, flag_ok)
% usage: 
%    [Veldata] = get_adcp_data(beg_time, end_time, index_file, data_path, vars);
%    returns data from the ps01 cruise for a specified time range into the 
%    Veldata struct.
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

Veldata = [];

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
   vars = []; % return all variables
end
if nargin<6
   flag_ok = 0;  % don't allow special overrides
end
if (end_time-beg_time)*24>48 & flag_ok~=1
   disp('Too much data requested, exceeds 48 hours');
   return
end

% check for existence of index file, if it exists load it
% otherwise return with an error and explaination.
if(exist(index_file) == 0)
   disp('Error:  Index file does not exist.\n');
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

% Find and sort files in range
yb = Index(iSet).yday_beg; ye = Index(iSet).yday_end;
ifil = find(yb<end_time & ye>beg_time);
if isempty(ifil)
   disp(['Index has no files for specified time range'])
   return
end
yb = yb(ifil); ye=ye(ifil);
fn=Index(iSet).filename; fn=fn(ifil);
if 0
    fzg = Index(iSet).zgrid(ifil,1:3);
    fz0 = fzg(1,1:3); % first z-grid -- all others must match this
    izn = find(fzg(:,1)~=fz0(1) | fzg(:,2)~=fz0(2) | fzg(:,3)~=fz0(3));
    if ~isempty(izn) 
        disp('Some ADCP files in range have unequal zgrids; Omitted.')
        fn(izn)=[]; yb(izn)=[]; ye(izn)=[];
    end
end
[x,iord] = sort(yb);
fyb = yb(iord); fye = ye(iord); fnam = fn(iord);

% value for	correcting profile-times to actual UTC (EG, sometimes off by 1 hr)
yday_ref = datenum(Set_params(iSet).year,1,1); % start of year
if isfield(Set_params(iSet), 'ping_seconds_offset')
    yd_off = Set_params(iSet).ping_seconds_offset / 86400; 
elseif isfield(Set_params(iSet), 'yday_seconds_offset')
    yd_off = Set_params(iSet).yday_seconds_offset / 86400; 
else
    yd_off = 0;
end

%% Read in ADCP files (gridded)
for ifil = 1:length(fyb)
   Vel = [];
   load( fullfile(data_path, fnam{ifil}) );
   if(isempty(Vel))
      error(['Error:loading file: ' data_path, fnam{ifil} ]);
   else
      Vel2(ifil) = Vel;
   end;
end;

if (length(vars) == 0)
    vars = {'ens_no','z_adcp','p_adcp','yday','bottomBT','u_wat','v_wat', 'w_wat', ...
            'u_shipBT','v_shipBT','latFP','lonFP','latLP','lonLP','ldayFP','ldayLP', ...
            'heading'};
    % vars = {'ens_no', 'z_adcp','p_adcp','pulselen','yday','bottomBT','u_wat','v_wat','w_wat', ...
    %       'u_shipBT','v_shipBT','latFP','lonFP','latLP','lonLP','ldayFP','ldayLP', ...
    %       'lPCtimeoff','lSMG','lDMG','heading','roll'};
end;

vars = vars(find(~strcmp(vars, 'z_adcp')));
vars = vars(find(~strcmp(vars, 'yday')));
len = length(vars);
vars = vars(find(~strcmp(vars, 'p_adcp')));
if( len ~= length(vars))
   p_flag = 1;
else 
   p_flag = 0;
end;

% initialize the variables in the new struct.
for i=1:length(vars)
    % first check that vars exist in the structures.
    if( isfield(Vel, vars{i}) )
        Veldata = setfield(Veldata, vars{i},[]);
    else
        disp(['Excluding invalid ADCP field=' vars{i}])
        vars{i} = '';
    end;
end;

Veldata.yday = [];
for j=1:length(Vel2)
    clear Vel
    Vel = Vel2(j);
    index = find(Vel.yday+yd_off >= beg_time & Vel.yday+yd_off <= end_time);
    if ~isempty(index)
        for i=1:length(vars)
            if ~isempty(vars{i})
                temp = ['[ Veldata.' vars{i} ', Vel.' vars{i} '(:,index)]'];
                temp2 = eval(temp);
                Veldata = setfield(Veldata, vars{i}, temp2);
            end
        end;
        Veldata.yday = [Veldata.yday  Vel.yday(index)+yd_off];
    end
end;
% Get depths (and pressures, if specified in vars-input) from first file
Veldata.z_adcp = Vel2(1).z_adcp';
if(p_flag == 1)
    Veldata.p_adcp = Vel2(1).p_adcp';
end;

