function [OVT1, OVT2] = Calc_ovts_grid(CSci, yd_pfb, yd_pfe)
% [OVT1, OVT2] = Calc_ovts_grid(CSci, yd_pfb, yd_pfe);
%   given 24-Hz CTD data (CSci) with yearday profile limits (yd_pfb,yd_pfe),
%   compute overturn statistics; Also return results on defined pressure grid.

%% filter 24-Hz data for sorting
SecBW = 2; % length of 4th-order Butterworth filter to lowpass pressure data
HzSamp = 24; % sampling frequency of pres,yday data
PtoZ = 1; % factor to convert pressures to depth in meters (2015: 1m=1db)
NumBW = ceil(SecBW*HzSamp/2)*2; 
NhfBW = NumBW/2; % half lengths
% Get filter coefs
if NumBW >= 2
    [bBW,aBW] = MHAButter(1/HzSamp ,SecBW);
end
%% find index limits, for filtering (profile -/+ 10 sec)
clear CTD
OVT1=[]; OVT2=[]; CTD=[]; ydS=[];
xtr = 10;
iFlt = find( CSci.yday_adj>=yd_pfb-(xtr/86400) & CSci.yday_adj<=yd_pfe+(xtr/86400) );
if length(iFlt) > NumBW*4 +(xtr*HzSamp) % ensure enough data
    % Filter profile variables
    CTD.T1=medfilt1(CSci.T1(iFlt),5); CTD.T2=medfilt1(CSci.T2(iFlt),5);
    CTD.S1=medfilt1(CSci.S1(iFlt),5); CTD.S2=medfilt1(CSci.S2(iFlt),5);
    CTD.Sg1=medfilt1(CSci.Sg1(iFlt),7); CTD.Sg2=medfilt1(CSci.Sg2(iFlt),7);
    % Filter pressure
    prS = CSci.Pr(iFlt);
    prLP = prS;
    ydS = CSci.yday_adj(iFlt);
    if NumBW >= 2 % Apply lowpass BW (2-way), if indicated
        prLP = filtfilt(bBW, aBW, prS);
    end
end
iPro = find(ydS>=yd_pfb & ydS<=yd_pfe); % Trim to profile limits
ix = find(iPro<=NhfBW | iPro>=length(ydS)-NhfBW); % exclude filter effects near ends
if ~isempty(ix)
    iPro(ix)=[];
end
if length(iPro) > NumBW*4 % enough data to proceed
    CTD.pr = prLP(iPro)' ; % dbar
    CTD.yday = ydS(iPro)';
    fn={'T1','T2','S1','S2','Sg1','Sg2'}; 
    for i=1:length(fn)
        eval(['CTD.' fn{i} ' = CTD.' fn{i} '(iPro)'';'])
    end
    % CTD.S1 = CTD.S1; CTD.S2 = CTD.S2; % PSU already, was *1000 for C.U.
    CTD.year=CSci.year;
    %% Now, Thorpe sort and compute overturn stats using densities
    clear SWOvt1 SWOvt2 PP
    %Setup PP, for Thorpe sort parameters
    PP=DefaultOverturnPP;
    %SWIMS Adjust 
    PP.cvar='';
    PP.pvar='pr';
    PP.drho=2e-4 / 2; % try for medfilt1'd T,S,Sigth
    PP.wh=1; % one at a time, using filtered 24-Hz data
    PP.threshold=0.1;
    PP.gridded='no';
    PP.plotit=0;
    PP.n = 5;
    % Do for both T,C pairs
    PP.tvar='T1'; PP.svar='S1';
    PP.cvar=''; PP.sgthvar='Sg1';
    [SWOvt1,OVT1]=FindOverturnsFCN2(CTD,PP);
    [SWOvt1,OVT1]=EpsFromOverturns2(SWOvt1,OVT1);
    %keyboard
    PP.tvar='T2'; PP.svar='S2';
    PP.cvar=''; PP.sgthvar='Sg2';
    [SWOvt2,OVT2]=FindOverturnsFCN2(CTD,PP);
    [SWOvt2,OVT2]=EpsFromOverturns2(SWOvt2,OVT2);
    clear SWOvt1 SWOvt2
end
clear CTD ydS prS prLP PP

