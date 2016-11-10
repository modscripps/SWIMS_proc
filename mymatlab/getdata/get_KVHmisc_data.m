function [KVHdata] = get_KVHmisc_data(beg_time, end_time, index_file, data_path, flag_ok)
% usage: 
%    [KVHdata] = get_KVHmisc_data(beg_time, end_time, index_file, data_path);
%    returns KVH data (from various serial sources) for a specified time range 
%    into the KVHdata struct.
%
% input:
%    beg_time   - lower limit of the time for which data is to be retrieved. (yeardays)
%    end_time   - upper limit of the time for which data is to be retrieved. (yeardays)
%    index_file - the file which indexes the data files.
%    data_path  - path to the data files which are indexed.
%
% output:
%    KVHdata    - a struct with fields of 'PCtimestamp', 'PCelapsetime', and other vectors
%       that depend on the data source.
% Dave W - 9/2002
KVHdata = [];

if nargin<3 | isempty(index_file)
	disp(['get_KVHmisc_data: Inputs beg_time,end_time,index_file are required.']);
    return
end
if beg_time<0 | end_time>370 | beg_time>=end_time
   disp(['get_KVHmisc_data: beg_time,end_time are out of range/order.']);
   return
end
if nargin<4
   data_path = []; % already in path (??)
end
if nargin<5
   flag_ok = 0;  % don't allow special overrides
end
if (end_time-beg_time)*24>48 & flag_ok~=1
   disp('get_KVHmisc_data: Too much data requested, exceeds 48 hours');
   return
end

% check for existence of index file, if it exists load it
% otherwise return with an error and explanation.
if(exist(index_file) == 0)
   disp('get_KVHmisc_data: Index file does not exist.\n');
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
   disp(['get_KVHmisc_data: Index does not cover specified time range'])
   return
end

% Find and sort files in range
yb = Index(iSet).yday_beg; ye = Index(iSet).yday_end;
ifil = find(yb<end_time & ye>beg_time);
if isempty(ifil)
   disp(['get_KVHmisc_data: Index has no files for specified time range'])
   return
end
yb = yb(ifil); ye=ye(ifil);
fn=Index(iSet).filename; fn=fn(ifil);
[x,iord] = sort(yb);
fyb = yb(iord); fye = ye(iord); fnam = fn(iord);

yday_ref = datenum(Set_params(iSet).year,1,1); % start of year

%% Read in KVH files 
for ifil = 1:length(fyb)
    KVH = []; LT = [];
    load( fullfile(data_path, fnam{ifil}) );
    if(isempty(KVH)) && ~isempty(LT)
        KVH = LT; clear LT % adapted to get LT data (feb-2008)
    end
    if isempty(KVH)
        error(['Error:loading file: ' data_path, fnam{ifil} ]);
    else
        Kvh2(ifil) = KVH;
    end;
end;

vars = fieldnames(KVH);

% initialize the variables in the new struct.
for i=1:length(vars)
   % first check that vars exist in the structures.
   if( isfield(KVH, vars{i}) )
       KVHdata = setfield(KVHdata, vars{i},[]);
   else
      disp('Error: File is missing expected field.');
      KVHdata = [];
      return
   end;
end;
% Gather data from files into output structure
for j=1:length(Kvh2)
    clear KVH
    KVH = Kvh2(j);
    index = find(KVH.PCtimestamp >= beg_time & KVH.PCtimestamp <= end_time);
    if ~isempty(index)
        for i=1:length(vars)
            temp = ['[ KVHdata.' vars{i} ', KVH.' vars{i} '(:,index)]'];
            temp2 = eval(temp);
            KVHdata = setfield(KVHdata, vars{i}, temp2);
        end;
    end
end;
