function [CTDsci] = get_SWIMS_SciData(beg_time, end_time, raw_index, raw_path, flag_ok, params)
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

%%%  ADD DATA CLEANUP!!! (e.g. c1/c2 comparison, etc.)

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
if nargin<5 | isempty(flag_ok)
    flag_ok = 0;  % don't allow special overrides
end
if nargin<6;
    params = []; % use defaults for lags, calibs, etc.
end
if (end_time-beg_time)*24>12 & flag_ok~=1
    warning('Too much data requested, exceeds 12 hours');
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
TLalfa = 0.03; TLbeta = 1/7; % Cond thermal mass
dox_tau0 = 0; dox_tauFt = 0; % Dox time constant
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
    %% params from feb05 tests:
    cond_lag=[.65]; TLbeta=[1/6.5]; TLalfa=[0.023];
    dox_lag=-2.0; dox_tau0=1.25; dox_tauFt=0;
    if yday > 71 % SWIMS 1
        % dox_lag=-3.3;cond_lag=-1.2; % ML04
        dox_lag = -0.8; cond_lag = -1.2;
        %% from feb05 tests:
        cond_lag=[-0.63]; TLbeta=[1/6.5]; TLalfa=[0.023];
        dox_lag=-2.5; dox_tau0=1.5;
        %Test:
        %dox_lag=-2.3;cond_lag=-0.7;
    end
else
    cond_lag = 0.75;
    dox_lag = -5; flu_lag = -2.5;
    pH_lag = -1.7;
end

% check for parameters passed in via func call
    pvars = {'cond_lag','dox_lag','flu_lag','pH_lag','h2s_lag', ...
        'TLalfa','TLbeta', 'dox_tau0','dox_tauFt'};
if isstruct(params)
    for i=1:length(pvars)
        if isfield(params,pvars{i})
            x = eval(['params.' pvars{i} ';']);
            if ~isempty(x)
                eval([pvars{i} ' = x;']);
            end
        end
    end
elseif strcmp(params,'nolags'),
    for i=1:length(pvars)
        eval([pvars{i} ' = 0;']);
    end
    disp('Applying no lags!');
end

% if time,thermal lags not specified for C2, use C1's:
if length(cond_lag)==1
    cond_lag = [cond_lag, cond_lag];
end
if length(TLalfa)==1
    TLalfa = [TLalfa, TLalfa];
end
if length(TLbeta)==1
    TLbeta = [TLbeta, TLbeta];
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
CTDsci.TLalfa = TLalfa;
CTDsci.TLbeta = TLbeta;
CTDsci.dox_lag = dox_lag;
CTDsci.dox_tau0 = dox_tau0;
CTDsci.dox_tauFt = dox_tauFt;
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
%Get cond (jan-2005: actual meas, no lags/corrections here)
[a,b,c,d,m]=read_Ccal_sbe_SWIMS(snc1,year,yday);
CTDsci.C1=C_sbe(SWraw.cfreq/256, a,b,c,d,m, CTDsci.T1, CTDsci.Pr);
[a,b,c,d,m]=read_Ccal_sbe_SWIMS(snc2,year,yday);
CTDsci.C2=C_sbe(SWraw.cfreq2/256, a,b,c,d,m, CTDsci.T2, CTDsci.Pr);
%% De-spiking of Temperature & Conductivity data.
% utilize sensor difference for outlier detection:
sC1=CTDsci.C1;sC2=CTDsci.C2;sT1=CTDsci.T1;sT2=CTDsci.T2;
[CTDsci.C1,CTDsci.C2,flag_c1,flag_c2] = two_sensor_cleanup(CTDsci.C1,CTDsci.C2);
[CTDsci.T1,CTDsci.T2,flag_t1,flag_t2] = two_sensor_cleanup(CTDsci.T1,CTDsci.T2);
keyboard

%%
clear S1 S2 T1 T2
%Compute salinity and density here. Check for 'Fresh' flags on C1,C2 first

if ~isempty(sparec1) & strcmpi(sparec1(1),'F')
    CTDsci.S1 = 0 * CTDsci.Pr;
    disp('Fresh water');
else
    clear CPs
    CPs.lag = cond_lag(1);
    CPs.alfa = TLalfa(1);
    CPs.beta = TLbeta(1);
%     S1 = salin_swims(CTDsci.C1, CTDsci.T1, CTDsci.Pr, CPs);
    % Conductivity time series need to be de-spiked at this point,
    % otherwise the spikes get distorted into strange things
    S1 = salin_swims_fft(CTDsci.C1, CTDsci.T1, CTDsci.Pr, CPs);
    CTDsci.S1 = S1.salin;
    clear S1
end
CTDsci.Th1 = sw_ptmp(CTDsci.S1, CTDsci.T1, CTDsci.Pr, 0);
CTDsci.Sg1 = sw_dens(CTDsci.S1, CTDsci.Th1, 0)-1000;
if ~isempty(sparec2) & strcmpi(sparec2(1),'F')
    CTDsci.S2 = 0 * CTDsci.Pr;
else
    clear CPs
    CPs.lag = cond_lag(2);
    CPs.alfa = TLalfa(2);
    CPs.beta = TLbeta(2);
%     S2 = salin_swims(CTDsci.C2, CTDsci.T2, CTDsci.Pr, CPs);
    S2 = salin_swims_fft(CTDsci.C2, CTDsci.T2, CTDsci.Pr, CPs);
    CTDsci.S2 = S2.salin;
    clear S2
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
    DoxCoefs.dox_lag = dox_lag;
    DoxCoefs.dox_tau0 = dox_tau0; DoxCoefs.dox_tauFt = dox_tauFt;
    % doxV_adj = shift( AtoD_SWIMS(SWraw.addata(ibN,:),SwimNo), dox_lag);
    doxV_adj = AtoD_SWIMS(SWraw.addata(ibN,:), SwimNo); % shift in Dox_sbe (jan05)
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



function [C1,C2,flag1,flag2] = two_sensor_cleanup(C1,C2);
% Perform spike detection and replacement using the two-sensor data
% (outliers are likely to affect each sensor independently!)

dC = C2 - C1;
dCm = mean(dC);
n = length(dC);
dCmax = sqrt(2*log(n))*std(dC); % maximum "natural" cond. difference between the sensors, based on the statistics
dCmax = 2*dCmax; % relax constraint a bit
% cf. Goring, D., V. Nikora (2002) Despiking Acoustic Doppler Velocimeter Data Journal of Hydraulic Engineering, Vol. 128, No. 1, p117- DOI: 10.1061/(ASCE)0733-9429(2002)128:1(117)
is = find(abs(dC-mean(dC))>dCmax); % outlier candidates
% which series is to blame? look at the second derivatives...
d2C1 = zeros(n,1);
d2C2 = zeros(n,1);
d2C1(2:n-1) = C1(3:n)-2*C1(2:n-1)+C1(1:n-2);
d2C2(2:n-1) = C2(3:n)-2*C2(2:n-1)+C2(1:n-2);
bad_in_C1 = abs(d2C1(is))>abs(d2C2(is));
is1 = is(bad_in_C1);
is2 = is(~bad_in_C1);
% add a couple of scans - deckbox shifting may have spread the spikes...
is1 = unique([is1;is1+1;is1-1]);
is1 = is1(is1>1 & is1<n);

is2 = unique([is2;is2+1;is2-1]);
is2 = is2(is2>1 & is2<n);

% fix C1 and C2
C1(is1) = C2(is1)-dCm; % use the other sensor to fill the gap
C2(is2) = C1(is2)+dCm; % use the other sensor to fill the gap
if nargout>2,
    flag1 = ones(n,1);
    flag2 = ones(n,1);
    flag1(is1) = 0;
    flag2(is2) = 0;
end