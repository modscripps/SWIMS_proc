function [ADdata] = get_ADCP_Any_vecs(beg_time, end_time, index_file, data_path, flag_ok)
%    [ADdata] = get_ADCP_Any_vecs(beg_time, end_time, index_file, data_path);
%    returns data from any matlab ADCP data in the SWIMS database.
%    That is, SWIMS mounted ADCPs raw ADUP,ADDN and processed VelUP,VelDN,
%    and vessel-mounted (via VMDas) ADCP and (e.g.) RDI_All .
%    get_ADCP_Any_vecs() EXCLUDES 2-dim arrays, to look at large numbers of
%         ensembles without depleting memory (e.g. for long mooring deployments)
% example:
%    AD = get_ADCP_Any(124.8, 124.9, ...
%       'C:\swims\ps02\indexes\ADDN_ps02_matfiles.mat', ...
%		'C:\swims\ps02\data_mat\ADDN\');
% input:
%    beg_time   - lower yearday limit for which data are to be retrieved.
%    end_time   - upper yearday limit for which data are to be retrieved.
%    index_file - the file which indexes the matlab data files (with '.mat').
%    data_path  - path to the matlab data files which are indexed.
%    (optional fifth argument = 1 to get more than 48 hrs of data)
% output:
%    ADdata - a structure with all fields common to ALL matlab data files
%      within the requested yearday limits.  If the depth grid changes within
%      this time period, all profile data are interpolated onto the grid used
%      by the first file (e.g. for BS03, vessel-based ADCP was switched
%      between the 150-kHz BroadBand and the 75-kHz Ocean Surveyor).
%    Fields always present are ADdata.yday, ADdata.ens_no and ADdata.z_adcp;
%      Others depend on the particular type of ADCP data.

% revision history:
%    20-sep-2003 - original version (get_ADupdn_data.m) for SWIMS-mounted ADCPs
%    27-may-2004 - generalized to work for VMDas-based data, and to
%           interpolate onto depth grid of earliest data in range - Dave W

ADdata = [];

if nargin<3 | isempty(index_file)
	disp(['get_ADCP_Any: Inputs beg_time,end_time,index_file are required.']);
    return
end
if beg_time<0 | end_time>370 | beg_time>=end_time
   disp(['get_ADCP_Any: beg_time,end_time are out of range/order.']);
   return
end
if nargin<4
   data_path = []; % already in path (??)
end
if nargin<5
   flag_ok = 0;  % don't allow special overrides
end
if (end_time-beg_time)*24>48 & flag_ok~=1
   disp('get_ADCP_Any: Too much data requested, exceeds 48 hours');
   return
end

% check for existence of index file, if it exists load it
% otherwise return with an error and explanation.
if(exist(index_file) == 0)
   disp('get_ADCP_Any:  Index file does not exist.');
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
   disp(['get_ADCP_Any: Index does not cover specified time range'])
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
   disp(['get_ADCP_Any: Index has no files for specified time range'])
   return
end

yb = yb(ifil); ye=ye(ifil);
fn=Index(iSet).filename; fn=fn(ifil);
[x,iord] = sort(yb);
fyb = yb(iord); fye = ye(iord); fnam = fn(iord);
fn = [];

yday_ref = datenum(Set_params(iSet).year,1,1); % start of year

%% Read in ADCP files 
clear AD, AD = [];
fil_one = 0; % set to 1 after first valid file is found
zgDIF = 0; % set to 1 if differences in depth grids exist
for ifil = 1:length(fyb)
    SN = []; XX = [];
    XX = load( fullfile(data_path, fnam{ifil}) );
    if ~isempty(XX), SN=fieldnames(XX); end
    if ~isempty(SN) & ~( strcmp(SN{1},'ADUP') | strcmp(SN{1},'ADDN') | ...
            strcmp(SN{1},'VelUP') | strcmp(SN{1},'VelDN') | ...
            strcmp(SN{1},'SWIMS_Vels') | strcmp(SN{1},'Vel') )
        disp(['Error: loading file: ' fullfile(data_path, fnam{ifil}) ]);
        return
    elseif ~isempty(SN)
        SN=SN{1};
        if strcmp(SN,'SWIMS_Vels') 
            XX.SWIMS_Vels.z_adcp = XX.SWIMS_Vels.depth;
            XX.SWIMS_Vels.yday = XX.SWIMS_Vels.yday_beg;
        end
        eval(['fn = fieldnames(XX.' SN ');']); % field names
        eval(['zg = XX.' SN '.z_adcp;']); % depth grid
        eval(['yddd = XX.' SN '.yday;']); % yeardays
        if length(zg)<2
            continue % for ensembles with no depth bins (oct-2004)
        end
        %keyboard
        if ~fil_one % first valid file
            eval(['AD = XX.' SN ';'])
            % Keep profile arrays from first file, for now
            fn0 = fn; % keep intersection of fields
            zg0 = zg; % if depth grids differ, use first file's
            fil_one = 1;
        else
            fn0c = fn0;
            zgEQ = isequal(zg, zg0); % detect unequal gridding
            % if only negligable differences, assume grids are same
            if ~zgEQ & size(zg)==size(zg0)
                dg = diff(zg(1:2)); dg0 =diff(zg0(1:2));
                zz = abs(zg(1)-zg0(1)) / dg;
                % check bin size and first bin depth, okay if close
                if abs(dg-dg0)/dg0 < 0.01 & zz < 0.02
                    zgEQ = 1; zg = zg0;
                end
            end
                
            for i=1:length(fn0)
                ok = eval(['isfield(XX.' SN ', ''' fn0{i} ''')']);
                if ok
                    eval(['FF = XX.' SN '.' fn0{i} ';'])
                    if ~strcmp(fn0{i}, 'z_adcp') & ~strcmp(fn0{i}, 'p_adcp') & ...
                            ~strcmp(fn0{i}, 'depth') & ~strcmp(fn0{i}, 'depSH') & ...
                            ~strcmp(fn0{i}, 'pulselen')
%                         if ~zgEQ & size(FF,2)==length(yddd) & ...
%                                 size(FF,1)==length(zg)
%                             x = FF;
%                             FF = interp1(zg,x, zg0); % interpolate onto first grid
%                             x = [];
%                             if ~zgDIF
%                                 zgDIF=1;
%                                 disp('WARNING: Unequal ADCP grids, interpolating profiles!')
%                             end
%                         end
                        % append here, but exclude profile arrays:
                        if ~(size(FF,2)==length(yddd) & size(FF,1)==length(zg))
                            eval(['AD.' fn0{i} ' = [AD.' fn0{i} ', FF ];'])
                        end
                        FF = [];
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
vars = vars(find( ~strcmp(vars, 'z_adcp') & ~strcmp(vars, 'p_adcp') & ...
    ~strcmp(vars, 'depth') & ~strcmp(vars, 'depSH') & ...
    ~strcmp(vars, 'pulselen') ));
ADdata.z_adcp = AD.z_adcp;
if isfield(AD,'p_adcp')
    ADdata.p_adcp = AD.p_adcp;
end
if isfield(AD,'depth')
    ADdata.depth = AD.depth;
end
if isfield(AD,'depSH')
    ADdata.depSH = AD.depSH;
end

index = find(AD.yday >= beg_time & AD.yday <= end_time);
if ~isempty(index)
    for i=1:length(vars)
        if ~isempty(vars{i})
            okvec = [' & size(AD.' vars{i} ', 1) < length(AD.z_adcp)']; % vectors only
            ok = eval(['size(AD.' vars{i} ', 2) == length(AD.yday)' okvec]);
            if ok  % right size to save and return
                temp = ['ADdata.' vars{i} ' = AD.' vars{i} '(:,index);'];
                eval(temp);
            end
        end
    end
end


