function CTDsci = calc_SWIMS_SciData(SWraw, params)
% usage: 
%  CTDsci = calc_SWIMS_SciData(SWraw, params);
%  Given 24-Hz raw data structure SWraw, applies calibrations to return
%   scientific data in structure CTDsci.
%   This routine was originally part of get_SWIMS_SciData.m, but was
%   spilt off in 2015 to allow alternate methods of generating SWraw.
%   Calibrations are applied to all active channels, creating fields:
%   Pr, T1,T2, C1,C2, Roll,Pitch, [SWIMS2: Dox,Flu,Obs,Obs2(oct/2011)];
%   and derived: S1,S2, Th1,Th2, Sg1,Sg2 (salin, theta, sigma_theta),
%   timestamps: yday_adj,SWIMStime,SWIMSfasttime, modCT,
%   scalars: cond_lag, year, SwimNo, params-fields.
%  NOTE: SWraw.Pr has already been computed, via Sav_rawCTs_SWIMS_wPR.m
%  params = structure of calibration-related parameters, if using values
%   other than the default ones specified below.
%   Dave W - 6/2015

CTDsci = [];

if nargin<2;
    params = []; % use defaults for lags, calibs, etc.
end

if ~isstruct(SWraw)
    warning('SWraw is not a raw data structure')
    return
end

year = SWraw.year;
SwimNo = SWraw.SWIMS_num;
% Use yearday at middle of data interval to determine sensors, calibrations
yday = (SWraw.yday_adj(1)+SWraw.yday_adj(end)) / 2;

%% for now, set lags by checking year,yearday
cond_lag=0; dox_lag=0; flu_lag=0; pH_lag=0; h2s_lag=0;
cond_tau=0; sal_fft=0; % Cond resp to match Temp, use time domain methods for salin
TLalfa = 0.03; TLbeta = 1/7; % Cond thermal mass
dox_tau0 = 0; dox_tauFt = 0; % Dox time constant
t1_lag=0; t2_lag=0; % Scans; to adjust C/T cross cabling for fresh water (pump activating)
if year<2002
    cond_lag = 1.350;
    %% from feb05 tests (swims1 w/dox, though):
    cond_lag=[-0.63]; TLbeta=[1/6.5]; TLalfa=[0.023];
elseif year==2002 && yday<180
    cond_lag = 1.350;
    dox_lag = -8.5; flu_lag = -2.5;
elseif year==2002 && yday>=180
    dox_lag = -6.2; flu_lag = -2.5;
elseif year==2003 && yday<180
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
elseif year==2007
    %% July 6, implemenent and use FFT method for salin
    sal_fft = [1, 0]; % use older method for secondary
    %% next three were for fft method, without first difn'ing
%     cond_lag = [1.78, 0.5]; % SWIMS3, per Andrey (then adjusted)
%     cond_tau = [0.033, 0]; 
%     TLbeta=[1/12.80, 1/6.5]; TLalfa=[0.015, 0.023];
    cond_lag = [1.934, 0.5]; % SWIMS3, per Andrey (for fft,fd'd)
    cond_tau = [0.0381, 0]; 
    TLbeta=[1/11.747, 1/6.5]; TLalfa=[0.0142, 0.023];
    %% for C2,T2 fft,fd'd: lag=1.989,tau=.0609,b=1/13.072,a=.0164
    dox_lag=-2.0; dox_tau0=1.25; dox_tauFt=0;
    flu_lag = -2.5;
elseif year==2009 % mc09
    sal_fft = [1, 1]; % use new method for both
    % Params per Andrey, run 05-mar-2010 pair#1, 09-feb-2010 pair#2 
    cond_lag = [1.924, 1.692]; % (for fft,fd'd)
    cond_tau = [0.0321, 0.0327]; 
    % Therm Lag, use previous for pair#1, feb-2010 for pair#2
    TLbeta=[1/11.747, 1/11.840]; TLalfa=[0.0142, 0.0196];
    dox_lag=-2.0; dox_tau0=1.25; dox_tauFt=0;
    flu_lag = -2.5;
elseif year==2011  || year==2013 % MORT (or wa_nliw_apr2013 for now)
    sal_fft = [1, 1]; % use new method for both
    cond_lag = [1.5,1.35]; % (for fft,fd'd)
    cond_tau = [0.0321, .03]; 
    % Therm Lag, use previous for pair#1, feb-2010 for pair#2
    TLbeta=[1/11.747, 1/11.840]; TLalfa=[0.0142, 0.0196];
    dox_lag=-2.0; dox_tau0=1.25; dox_tauFt=0;
    flu_lag = -2.5 -1;
elseif year == 2016
    cond_lag = 0.75;
    dox_lag = -5;
    flu_lag = -2.5;
    pH_lag = -1.7;
else
    cond_lag = 0.75;
    dox_lag = -5; flu_lag = -2.5;
    pH_lag = -1.7;
end

% check for parameters passed in via func call
if isstruct(params)
    pvars = {'cond_lag','dox_lag','flu_lag','pH_lag','h2s_lag', ...
        'cond_tau','sal_fft','TLalfa','TLbeta', ...
        'dox_tau0','dox_tauFt'};
    for i=1:length(pvars)
        if isfield(params,pvars{i})
            x = params.(pvars{i});
            if ~isempty(x)
                eval([pvars{i} ' = x;']);
            end
        end
    end
end

% if time,thermal lags not specified for C2, use C1's:
if length(cond_lag)==1
    cond_lag = [cond_lag, cond_lag];
end
if length(cond_tau)==1
    cond_tau = [cond_tau, cond_tau];
end
if length(sal_fft)==1
    sal_fft = [sal_fft, sal_fft];
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
[sndox,sparedox]=GetSWIMSConfig('dox',year,yday); % 'sparedox' flags cal type (2011)

if SwimNo==2
    [snflu,spare]=GetSWIMSConfig('flu',year,yday);
    [snobs,spare]=GetSWIMSConfig('obs',year,yday);
    [snpH,spare]=GetSWIMSConfig('pH',year,yday);
    [snH2S,spare]=GetSWIMSConfig('H2S',year,yday);
    [snobs2,spare]=GetSWIMSConfig('obs2',year,yday); % OBS-3+ 'Low (NTU) Signal', oct/2011
end

% Save timestamps, other info
CTDsci.yday_adj = SWraw.yday_adj; % uses elapsed seconds + origin
CTDsci.SWIMStime = SWraw.SWIMStime; % based on integer seconds
CTDsci.SWIMSfasttime = SWraw.SWIMSfasttime; % elapsed seconds, float (msec)
CTDsci.modCT = SWraw.modCT;
CTDsci.Pr = SWraw.Pr; % Already done (calib'd) along with raw counts
CTDsci.status = SWraw.status;
% calib params
CTDsci.cond_lag = cond_lag;
CTDsci.cond_tau = cond_tau;
CTDsci.sal_fft = sal_fft;
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


% Check flags for switched C,T cabling (e.g., to activate pumps in fresh water)
if strcmpi(sparet1, 'c1')
    tmp = SWraw.tfreq; SWraw.tfreq = SWraw.cfreq; SWraw.cfreq = tmp;
    if year<2007
        t1_lag=2
    else % set deckbox to no delay
        t1_lag = 0
    end
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
%%% NOTE that as of jul 2007, deckbox shifts set=0 for SWIMS3 %%%
CTDsci.T1 = shift(T1,t1_lag); 
CTDsci.T2 = shift(T2,t2_lag);
%Get cond (jan-2005: actual meas, no lags/corrections here)
[a,b,c,d,m]=read_Ccal_sbe_SWIMS(snc1,year,yday); Copt='abcdm';
if isnan(m)
    Copt='ghij'; % coeffs for SBE's 'recommended' formula
end
CTDsci.C1=C_sbe(SWraw.cfreq/256, a,b,c,d,m, CTDsci.T1, CTDsci.Pr, Copt);
[a,b,c,d,m]=read_Ccal_sbe_SWIMS(snc2,year,yday); Copt='abcdm';
if isnan(m)
    Copt='ghij'; % coeffs for SBE's 'recommended' formula
end
CTDsci.C2=C_sbe(SWraw.cfreq2/256, a,b,c,d,m, CTDsci.T2, CTDsci.Pr, Copt);

% For now, remove entire scan with any bad C,T,Pr data
ie = find(isnan(CTDsci.T1+CTDsci.T2+CTDsci.C1+CTDsci.C2+CTDsci.Pr));
if ~isempty(ie)
    disp(['Removing ' int2str(length(ie)) ' bad (CTP) scans'])
    SWraw.addata(:,ie) = [];
    CTDsci.T1(ie) = [];
    CTDsci.T2(ie) = [];
    CTDsci.C1(ie) = [];
    CTDsci.C2(ie) = [];
    CTDsci.Pr(ie) = [];
    CTDsci.yday_adj(ie) = [];
    CTDsci.SWIMStime(ie) = [];
    CTDsci.SWIMSfasttime(ie) = [];
    CTDsci.modCT(ie) = [];
    CTDsci.status(ie) = [];
end
% Later ... clean spikes in C,T:
% %% De-spiking of Temperature & Conductivity data.
% % utilize sensor difference for outlier detection:
% sC1=CTDsci.C1;sC2=CTDsci.C2;sT1=CTDsci.T1;sT2=CTDsci.T2;
% [CTDsci.C1,CTDsci.C2,flag_c1,flag_c2] = two_sensor_cleanup(CTDsci.C1,CTDsci.C2);
% [CTDsci.T1,CTDsci.T2,flag_t1,flag_t2] = two_sensor_cleanup(CTDsci.T1,CTDsci.T2);
% keyboard

clear S1 S2 T1 T2
%Compute salinity and density here. Check for 'Fresh' flags on C1,C2 first
%%
if ~isempty(sparec1) && strcmpi(sparec1(1),'F')
    CTDsci.S1 = 0 * CTDsci.Pr;
else
    clear CPs
    CPs.lag = cond_lag(1);
    CPs.tau = cond_tau(1);
    CPs.alfa = TLalfa(1);
    CPs.beta = TLbeta(1);
    if sal_fft(1)
        S1 = salin_swims_fft(CTDsci.C1, CTDsci.T1, CTDsci.Pr, CPs);
    else
        S1 = salin_swims(CTDsci.C1, CTDsci.T1, CTDsci.Pr, CPs);
    end
    CTDsci.S1 = S1.salin;
%     clear S1    
end
CTDsci.Th1 = sw_ptmp(CTDsci.S1, CTDsci.T1, CTDsci.Pr, 0);
CTDsci.Sg1 = sw_dens(CTDsci.S1, CTDsci.Th1, 0)-1000;
if ~isempty(sparec2) && strcmpi(sparec2(1),'F')
    CTDsci.S2 = 0 * CTDsci.Pr;
else
    clear CPs
    CPs.lag = cond_lag(2);
    CPs.tau = cond_tau(2);
    CPs.alfa = TLalfa(2);
    CPs.beta = TLbeta(2);
    if sal_fft(2)
        S2 = salin_swims_fft(CTDsci.C2, CTDsci.T2, CTDsci.Pr, CPs);
    else
        S2 = salin_swims(CTDsci.C2, CTDsci.T2, CTDsci.Pr, CPs);
    end
    CTDsci.S2 = S2.salin;
%     clear S2  
end
CTDsci.Th2 = sw_ptmp(CTDsci.S2, CTDsci.T2, CTDsci.Pr, 0);
CTDsci.Sg2 = sw_dens(CTDsci.S2, CTDsci.Th2, 0)-1000;

% A-D data
RPver = SwimNo;
if SwimNo==2 && year>2004 % SWIMS 3, deck test 21-Feb-2008
    RPver = 3; % roll sense is opposite SWIMS2, pitch is same
end
CTDsci.Roll = roll_Swims( AtoD_SWIMS(SWraw.addata(1,:),SwimNo), RPver);
CTDsci.Pitch = pitch_Swims( AtoD_SWIMS(SWraw.addata(2,:),SwimNo), RPver);
if ~isnan(snalt)
    ALTcoefs = read_AltimCal_SWIMS(snalt,year,yday);
    CTDsci.Alt = Altim_Swims( AtoD_SWIMS(SWraw.addata(3,:),SwimNo), ...
        ALTcoefs.Mo, ALTcoefs.Bo);
end 
% Oxygen, special for SWIMS 1, ML04
if ~isnan(sndox)
    ibN = 4;
    if SwimNo==1 && isnan(snalt)
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
        if year==2004 && yday>59 && yday<65.8333 % switched sensor and to 30X to fix
            ix = find(CTDsci.Flu<.06);
            CTDsci.Flu(ix) = 10 * CTDsci.Flu(ix);
        end
    end
    if ~isnan(snobs)
        OBScoefs = read_ObsCal_SWIMS(snobs,year,yday);
        CTDsci.Obs = Obs_Swims( AtoD_SWIMS(SWraw.addata(7,:),SwimNo), ...
            OBScoefs.Mo, OBScoefs.Bo, OBScoefs.Qo);
    end
    if ~isnan(snobs2)
        OBS2coefs = read_ObsCal_SWIMS(snobs2,year,yday);
        CTDsci.Obs2 = Obs_Swims( AtoD_SWIMS(SWraw.addata(8,:),SwimNo), ...
            OBS2coefs.Mo, OBS2coefs.Bo, OBS2coefs.Qo);
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

