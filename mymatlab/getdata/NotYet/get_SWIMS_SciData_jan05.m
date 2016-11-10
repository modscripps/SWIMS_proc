function [CTDsci] = get_SWIMS_SciData(beg_time, end_time, raw_index, raw_path, flag_ok)
% usage: 
%  [CTDsci] = get_SWIMS_SciData(beg_time, end_time, index_file, raw_path);
%  Returns 24-Hz data for a specified time (yearday) range into the 
%   CTDsci structure.  Calibrations are applied to all active channels:
%   Structure fields are: Pr, T1,T2, C1,C2, Roll,Pitch, [SWIMS2: Dox,Flu,Obs];
%   derived: S1,S2, Th1,Th2, Sg1,Sg2 (salin, theta, sigma_theta),
%   timestamps: yday_adj,SWIMStime,SWIMSfasttime, modCT,
%   scalars: cond_lag, year, SwimNo.
% Inputs are yearday range, Raw-Swims data index file, Path to Raw data;
%
% EG savepath='c:/swims/ps02'; data_dir = fullfile(savepath,'data_mat');
%    IndFld = fullfile(savepath, 'indexes');
%    SWraw = get_SWIMS_RawData(118.7, 118.9, ...
%        fullfile(IndFld, 'CTD_ps02_matfiles.mat'), ...
%        fullfile(data_dir, 'CTD') );
% DPW - apr/2002

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

%% for now, set lags by checking year,yearday
cond_lag=0; dox_lag=0; flu_lag=0; pH_lag=0; h2s_lag=0;
t1_lag=0; t2_lag=0; % Scans; to adjust C/T cross cabling for fresh water (pump activating)
if year<2002
    cond_lag = 1.350;
elseif year==2002 & yday<180
    cond_lag = 1.350;
    dox_lag = -8.5; flu_lag = -2.5;
elseif year==2002 & yday>=180
    dox_lag = -6.2; flu_lag = -2.5;
elseif year==2003 & yday<180
    dox_lag = -6.2; flu_lag = -2.5;
    pH_lag = -1.7; h2s_lag = -11;
elseif year==2004
    dox_lag = -2.5; flu_lag = -2.3; 
    cond_lag = 0.2;
    pH_lag = -1.3;
    if yday > 71 % SWIMS 1
        % dox_lag=-3.3;cond_lag=-1.2; % ML04
        dox_lag = -0.8; cond_lag = -1.2;
        %Test:
        %dox_lag=-2.3;cond_lag=-0.7;
    end
else
    cond_lag = 0.75;
    dox_lag = -5; flu_lag = -2.5;
    pH_lag = -1.7;
end

dox_lag = dox_lag * 24; % convert to scans
flu_lag = flu_lag * 24; % convert to scans
pH_lag = pH_lag * 24; % convert to scans
h2s_lag = h2s_lag * 24; % convert to scans

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
CTDsci.status = SWraw.status;

% Check flags for switched C,T cabling (e.g., to activate pumps in fresh water)
if strcmpi(sparet1, 'c1')
    tmp = SWraw.tfreq; SWraw.tfreq = SWraw.cfreq; SWraw.cfreq = tmp;
    t1_lag=2
end
if strcmpi(sparet2, 'c2')
    tmp = SWraw.tfreq2; SWraw.tfreq2 = SWraw.cfreq2; SWraw.cfreq2 = tmp;
    t2_lag=2
end
clear tmp
%Get temperature 
[a,b,c,d,fo]=read_Tcal_sbe_SWIMS(snt1,year,yday);
T1=T_sbe(SWraw.tfreq/256, a,b,c,d,fo);
[a,b,c,d,fo]=read_Tcal_sbe_SWIMS(snt2,year,yday);
T2=T_sbe(SWraw.tfreq2/256, a,b,c,d,fo);
% Shift temperature, in case of cross-cabling with cond for fresh-water pumping,
%  to adjust for .073-s shift by deck-unit (SWIMS2)
CTDsci.T1 = shift(T1,t1_lag); 
CTDsci.T2 = shift(T2,t2_lag);
%Get cond 
[a,b,c,d,m]=read_Ccal_sbe_SWIMS(snc1,year,yday);
C1=C_sbe(SWraw.cfreq/256, a,b,c,d,m, CTDsci.T1, CTDsci.Pr);
[a,b,c,d,m]=read_Ccal_sbe_SWIMS(snc2,year,yday);
C2=C_sbe(SWraw.cfreq2/256, a,b,c,d,m, CTDsci.T2, CTDsci.Pr);
%Shift cond for lineup with temperature
%keyboard

CTDsci.C1=shift(C1,cond_lag);
CTDsci.C2=shift(C2,cond_lag);
clear S1 S2 T1 T2
%Compute salinity and density here. Check for 'Fresh' flags on C1,C2 first
%%
if ~isempty(sparec1) & strcmpi(sparec1(1),'F')
    CTDsci.S1 = 0 * CTDsci.Pr;
else
    CTDsci.S1 = salinityfcn(CTDsci.C1, CTDsci.T1, CTDsci.Pr);
end
CTDsci.Th1 = sw_ptmp(CTDsci.S1, CTDsci.T1, CTDsci.Pr, 0);
CTDsci.Sg1 = sw_dens(CTDsci.S1, CTDsci.Th1, 0)-1000;
if ~isempty(sparec2) & strcmpi(sparec2(1),'F')
    CTDsci.S2 = 0 * CTDsci.Pr;
else
    CTDsci.S2 = salinityfcn(CTDsci.C2, CTDsci.T2, CTDsci.Pr);
end
CTDsci.Th2 = sw_ptmp(CTDsci.S2, CTDsci.T2, CTDsci.Pr, 0);
CTDsci.Sg2 = sw_dens(CTDsci.S2, CTDsci.Th2, 0)-1000;

% A-D data
CTDsci.Roll = roll_Swims( AtoD_SWIMS(SWraw.addata(1,:),SwimNo), SwimNo); %, 1);
CTDsci.Pitch = pitch_Swims( AtoD_SWIMS(SWraw.addata(2,:),SwimNo), SwimNo); %, 1);
if ~isnan(snalt)
    ALTcoefs = read_AltimCal_SWIMS(snalt,year,yday);
    CTDsci.Alt = Altim_Swims( AtoD_SWIMS(SWraw.addata(3,:),SwimNo), ...
        ALTcoefs.Mo, ALTcoefs.Bo);
end 
% Oxygen, special for SWIMS 1, ML04
if ~isnan(sndox)
    ibN = 4;
    if SwimNo==1 & isnan(snalt)
        ibN = 3;
    end
    DoxCoefs = read_DoxCal_sbe_SWIMS(sndox,year,yday);
    doxV_adj = shift( AtoD_SWIMS(SWraw.addata(ibN,:),SwimNo), dox_lag);
    CTDsci.Dox = Dox_sbe( doxV_adj, CTDsci.T1, CTDsci.S1, CTDsci.Pr, DoxCoefs);
end

if SwimNo==2
    if ~isnan(snflu)
        FLUcoefs = read_FluCal_SWIMS(snflu,year,yday);
        fluV_adj = shift( AtoD_SWIMS(SWraw.addata(5,:),SwimNo), flu_lag);
        CTDsci.Flu = Fluoro_Swims( fluV_adj, FLUcoefs.Mo, FLUcoefs.Bo);
        % Kludge for ML04 behavior, jumping from 10X (cable) to 1X gain
        if year==2004 & yday>59 & yday<65.8333 % switched sensor and to 30X to fix
            ix = find(CTDsci.Flu<.06);
            CTDsci.Flu(ix) = 10 * CTDsci.Flu(ix);
        end
    end
    if ~isnan(snobs)
        OBScoefs = read_ObsCal_SWIMS(snobs,year,yday);
        CTDsci.Obs = Obs_Swims( AtoD_SWIMS(SWraw.addata(7,:),SwimNo), ...
            OBScoefs.Mo, OBScoefs.Bo);
    end
    if ~isnan(snpH)
        pHCoefs = read_pHCal_SWIMS(snpH,year,yday);
        pHV_adj = shift( AtoD_SWIMS(SWraw.addata(8,:),SwimNo), pH_lag);
        CTDsci.pH = pH_Swims( pHV_adj, CTDsci.T1, pHCoefs);
    end
    if ~isnan(snH2S)
        H2SCoefs = read_H2SCal_SWIMS(snH2S,year,yday);
        H2SV_adj = shift( AtoD_SWIMS(SWraw.addata(6,:),SwimNo), h2s_lag);
        %% check this,  first output is total, 2nd is H2S
        [CTDsci.Sulf, CTDsci.H2S] = ...
            H2S_Swims( H2SV_adj, CTDsci.T1, CTDsci.pH, CTDsci.Pr, H2SCoefs);
    end
end

% Subsample data to fall in requested yearday range 
vars = {'Pr','T1','T2','C1','C2','S1','S2','Th1','Th2','Sg1','Sg2', ...
        'Roll','Pitch','Dox','Flu','Obs','pH','Sulf','H2S','Alt','status'};
for i=1:length(vars)
   % make sure vars{i} exists for this data.
   if ( isfield(CTDsci, vars{i}) )
       eval(['CTDsci.' vars{i} ' = CTDsci.' vars{i} '(inRng);'])
   end
end
