function [CTDdata] = get_swims_data_combo(beg_time, end_time, index_file, data_path, vars, flag_ok)
% usage: 
%    [CTDdata] = get_swims_data(beg_time, end_time, index_file, data_path, vars);
%    returns data from the ps01 cruise for a specified time range into the 
%    CTDdata struct.
%
% example:
%    stuff = get_swims_data(98, 99, '/a/esp/swims_data/swims/ps01/griddata/', ...
%			'/home/zach/work/ps01/SWIMS_Log_ps01.mat', {'s1', 'th1'})
% 
% input:
%    beg_time   - lower limit of the time for which data is to be retrieved. (yeardays=prof_start)
%    end_time   - upper limit of the time for which data is to be retrieved. (yeardays=prof_start)
%    index_file - the file which indexes the data files. created by make_ps01_swims_index.m
%    data_path  - path to the data files which are indexed.
%    vars       - list of the vars of interest on the data files. ex: { 's1', 'th1'}
%                 if {} is passed in, the function returns all of the variables.
%
% output:
%    CTDdata    - a struct with fields of the requested variables plus 'z_ctd' and 'yday' 
%                 which are always returned.

% revision history:
%    8.17.01    - original version.              zach frazier
%    1.22.02    - Use revised file index, other stuff - Dave W

CTDdata = [];

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
   vars = []; % return all variables
end
if nargin<6 | isempty(flag_ok)
   flag_ok = 0;  % don't allow special overrides
end
if (end_time-beg_time)*24>48 & flag_ok<1
   warning('Too much data requested, exceeds 48 hours');
end
%keyboard
% check for existence of index file, if it exists load it
% otherwise return with an error and explaination.
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
   error(['Index does not cover specified time range'])
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
fzg = Index(iSet).zgrid(ifil,1:3);
% If requested, form union of differing grids (with coarsest spacing)
ReGrid = 0;
zg_min = min(fzg(:,1));
zg_del = max(fzg(:,2));
zg_max = max(fzg(:,3));
ZG_all = [zg_min:zg_del:zg_max];
fz0 = fzg(1,1:3); % first z-grid -- default = all others must match this
izn = find(fzg(:,1)~=fz0(1) | fzg(:,2)~=fz0(2) | fzg(:,3)~=fz0(3));
if ~isempty(izn)
    if flag_ok<2
        disp('Some SWIMS files in range have unequal z-grids; Omitted.')
        fn(izn)=[]; yb(izn)=[]; ye(izn)=[];
    else
        ReGrid = 1; % interpolate all onto ZG_all
    end
end
[x,iord] = sort(yb);
fyb = yb(iord); fye = ye(iord); fnam = fn(iord);

% value for	correcting profile-times to actual UTC (EG, sometimes off by 1 hr)
yday_ref = datenum(Set_params(iSet).year,1,1); % start of year
yd_off = Set_params(iSet).yday_seconds_offset / 86400; 

%% Read in SWIMS files (gridded), correct for midnight problem
for ifil = 1:length(fyb)
    SWIMSgrid = [];
    load( fullfile(data_path, fnam{ifil}) );
    if(isempty(SWIMSgrid))
        error(['Error:loading file: ' data_path, fnam{ifil} ]);
    else
        ix = find(SWIMSgrid.yday<SWIMSgrid.yday(1));
        SWIMSgrid.yday(ix) = SWIMSgrid.yday(ix)+1;
        fn = fieldnames(SWIMSgrid);
        %keyboard
        for i=1:length(fn)
            eval(['Ctd2(ifil).' fn{i} ' = SWIMSgrid.' fn{i} ';']);
        end
        %Ctd2(ifil) = SWIMSgrid;
    end;
end;

if length(vars) == 0
    varDEF = 1; % default: return all fields (profile attributes and gridded data)
    vars = fieldnames(SWIMSgrid); % use fields from last grid-file in yearday range
else
    varDEF = 0; % specified variables requested by user
end
vars = vars(find(~strcmp(vars, 'z')));
vars = vars(find(~strcmp(vars, 'p')));
vars = vars(find(~strcmp(vars, 'yday')));

% initialize the variables in the new struct.
for i=1:length(vars)
   % first check that vars exist in the structures.
   if( isfield(SWIMSgrid, vars{i}) )
       CTDdata = setfield(CTDdata, vars{i},[]);
   else
       if ~varDEF
           disp(['Excluding invalid SWIMS grid field=' vars{i}])
       end
       vars{i} = '';
   end
end

% ProVars = []; % Find fields that are profile quantities (for regridding)
% lenZ = length(SWIMSgrid.z);
% for i=1:length(vars)
%     lenV = eval(['size(SWIMSgrid.' vars{i} ',1)']);
%     if lenV==lenZ & isempty(find( strcmp(vars{i}, {'z','p'}) ))
%         % add to required fields, and fill with NaNs
%         ProVars{end+1} = vars{i};
%     end
% end

CTDdata.yday = [];
for j=1:length(Ctd2)
    clear CTD
    CTD = Ctd2(j);
    index = find(CTD.yday+yd_off >= beg_time & CTD.yday+yd_off <= end_time);
    if ~isempty(index)
        for i=1:length(vars)
            if ~isempty(vars{i})
                temp = ['[ CTDdata.' vars{i} ', CTD.' vars{i} '(:,index)]'];
                temp2 = eval(temp);
                CTDdata = setfield(CTDdata, vars{i}, temp2);
            end
        end;
        CTDdata.yday = [CTDdata.yday  CTD.yday(index)+yd_off];
    end
end;
% Get depths from first file
CTDdata.z_ctd = Ctd2(1).z';
CTDdata.p_ctd = CTDdata.z_ctd/100;

