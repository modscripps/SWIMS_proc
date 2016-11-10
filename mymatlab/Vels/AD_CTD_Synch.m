% AD_CTD_Synch.m  -  Find clock offsets between SWIMS ADCPs and CTD.
%  Mainly based on offset to ADDN water pings, with constant offsets
%  between ADDN water and BT pings, and between ADDN and ADUP water pings.
% Need to set variables: cruise, yd_b, yd_e  before running

% called with args from Run_SWIMS_sync.m
% cruise = 'philex08';
% yd_b = 298.3; yd_e = 298.32; % hc03 test
set_swims_paths

ex_yday = 20/86400; % extra data at ends

% Ping sequence and timing varies with cruise, BTing or not.
% Based on lab tests before home02, most reliable setting was
% to have ADUP as Master, ADDN as Slave:  This yielded sequence
% DnBT, DnW, UpW, with predictable UpW-DnBT timing, but with
% minor DnW-DnBT and UpW-DnW variations O(0.05-s) depending on
% BT conditions (which are not reproduceable in the lab).
%   The settings used before these tests (of Jul 2002) had
% ADDN as Master, ADUP as Slave, with seq = UpW, DnBT, DnW,
% reliable UpW-DnBT, but with DnW-DnBT and UpW-DnW
% variations Order(0.30-s) in the lab setting.
%   For periods with BT turned off, DnW, UpW; DnW, UpW; ...
% repeated every 1.05-s during bs03, which was sometimes too fast
% for data acquisition to keep up; Will try 2x0.6s for ML04.
% DelDNtoBT = (yday_BT-yday_DN), DelDNtoUP = (yday_UP-yday_DN) .
switch lower(cruise)
    case 'ps02' % sequence UpW, DnBT, DnW
        DelDNtoBT = -0.25;
        DelDNtoUP = -0.67;
    case 'home02'
        DelDNtoBT = -0.43;
        DelDNtoUP = 0.40; % SW = 08000
    case 'bs03'
        if yd_b>76.4 & yd_b<82.0 % period with no BTing
            DelDNtoBT = -0.1; % doesn't matter, not used
            DelDNtoUP = 0.51;
        else % while BTing
            DelDNtoBT = -0.43;
            DelDNtoUP = 0.40;
        end
    case 'hc03'
        DelDNtoBT = -0.44;
        DelDNtoUP = 0.47; % SW = 09000
    case 'ml04'
        if yd_b>59 & yd_e<91
            DelDNtoBT = -0.1; % doesn't matter, not used
            DelDNtoUP = 0.61; % With BT off in deep water
        else
            DelDNtoBT = -0.44;
            DelDNtoUP = 0.47; % SW = 09000
        end
    case 'stf07'  % now, SWIMS3
        DelDNtoBT = -0.1; % doesn't matter, not used
        DelDNtoUP = 0.61; % With BT off in deep water
    case 'philex08'
        if yd_e<44 % BT off
            DelDNtoBT = -0.1; % doesn't matter, not used
            DelDNtoUP = 0.61; % With BT off in deep water
        elseif yd_b>=44 && yd_e<999 % BT on
            DelDNtoBT = -0.44;
            DelDNtoUP = 0.47; % SW = 09000
        end
    case 'mc09'
        DelDNtoBT = -0.44;
        DelDNtoUP = 0.47; % SW = 09000
    case 'mort'
        DelDNtoBT = -0.44;
        DelDNtoUP = 0.47; % SW = 09000
    otherwise
        DelDNtoBT = -0.44;
        DelDNtoUP = 0.47; % SW = 09000
        % for no BTing after 2003:
        % DelDNtoBT = -0.1; % doesn't matter, not used
        % DelDNtoUP = 0.61;
end
DelDNtoBT = DelDNtoBT / 86400; % yeardays
DelDNtoUP = DelDNtoUP / 86400; % yeardays

yd_DN = []; yd_UP = [];
EnsDNtoUP = 0; % offset between ADUP and ADDN ensembles (UP-DN)
% Get needed fields from Downward ADCP
Vdn = get_ADupdn_data(yd_b-ex_yday, yd_e+ex_yday, ...
    ['/Users/ecfine/Documents/MATLAB/swims/ArcticMix/indexes/VelDN_' cruise '_matfiles.mat'], ...
    ['/Users/ecfine/Documents/MATLAB/swims/' cruise '/data_mat/VelDN']);
if ~isempty(Vdn)
    W_DN = median(Vdn.w(2:5,:));
    yd_DN = Vdn.yday;
    ens_DN = Vdn.ens_no;
    tip_DN = 1.4*Vdn.roll;
    W_BT = Vdn.wBT;
    yd_BT = yd_DN + DelDNtoBT;
    if isfield(Vdn, 'depth_xducer')
        z_DN = Vdn.depth_xducer;
    else
        z_DN = NaN*yd_DN;
    end
    clear Vdn; Vdn = 1; % for later logic
else
    disp('WARNING: No DN data')
    EnsDNtoUP = NaN;
end

% Get needed fields from Upward ADCP
Vup = get_ADupdn_data(yd_b-ex_yday, yd_e+ex_yday, ...
    ['/Users/ecfine/Documents/MATLAB/swims/' cruise '/indexes/VelUP_' cruise '_matfiles.mat'], ...
    ['/Users/ecfine/Documents/MATLAB/swims/' cruise '/data_mat/VelUP']);
if ~isempty(Vup)
    W_UP = median(Vup.w(2:5,:));
    yd_UP = Vup.yday;
    ens_UP = Vup.ens_no;
    tip_UP = 1.4*Vup.pitch; % should = DN roll = tip_DN
    if isfield(Vup, 'depth_xducer')
        z_UP = Vup.depth_xducer;
    else
        z_UP = NaN*yd_UP;
    end
    if isempty(Vdn)
        W_DN=NaN*W_UP; yd_DN=yd_UP; ens_DN=ens_UP; tip_DN=NaN*tip_UP;
        W_BT=W_DN; yd_BT=yd_DN; z_DN=z_UP;
    end
else 
    disp('WARNING: No UP data')
    W_UP=NaN*W_DN; yd_UP=yd_DN; ens_UP=ens_DN; tip_UP=NaN*tip_DN; z_DN=z_UP;
end
clear Vup

if isempty([yd_DN yd_UP])
    disp('No ADDN or ADUP data for this yearday range, quitting!')
    return
end

% Verify ensemble number synch between UP and DN
EnsDNtoUP = 0; % offset between ADUP and ADDN ensembles (UP-DN)
Fens = figure;
Fsyn = figure;
Fzad = 0; 
if ~isempty(~isnan(z_UP)) && ~isempty(~isnan(z_DN))
    Fzad = figure; % also check by comparing ADCP pr's (SWIMS3)
end
Fxtra=figure;
EnCont = 1;
if isnan(EnsDNtoUP)
    EnCont = 0; % skip this part if only UP
end
ens_DN_sav = ens_DN;
while EnCont
    figure(Fens), clf
    ens_UP_adj = ens_UP - EnsDNtoUP;
    % Here, make sure ensemble no.s are monotonic
    xeD = [find(diff(ens_DN)<1) length(ens_DN)];
    xeU = [find(diff(ens_UP_adj)<1) length(ens_UP_adj)];
    mxeD = NaN;
    for i=1:length(xeD)-1
        mxeD(i) = ens_DN(xeD(i)); % keep these to use for ADUP to keep ens_nos synched
        % force to be monotonic
        ens_DN([xeD(i)+1:xeD(i+1)]) = ens_DN([xeD(i)+1:xeD(i+1)]) + mxeD(i)+50;
    end
    for i=1:min(length(xeU),length(xeD))-1
        % force to be monotonic (use ADDN 'offsets' to keep synched)
        ens_UP_adj([xeU(i)+1:xeU(i+1)]) = ens_UP_adj([xeU(i)+1:xeU(i+1)]) + mxeD(i)+50;
    end
    if length(xeU) > length(xeD)
        ens_UP_adj([xeU(length(xeD))+1:end]) = inf; % beyond last ADDN ensemble
    end
    en0 = ens_DN(1);
    plot(ens_DN-en0,tip_DN,'r-',ens_DN-en0,tip_DN,'r.', ...
        ens_UP_adj-en0,tip_UP,'g:',ens_UP_adj-en0,tip_UP,'g.'), hold on, zoom on,grid on
    title(['UP(gn)--DN(rd) ens\_no = ' int2str(EnsDNtoUP) ', begin yd=' num2str(yd_DN(1))])
    xlabel(['Ensemble Number (DN) -- ' int2str(en0)])
    ylabel(['ADCPs Tip Angle'])

    if Fzad
        figure(Fzad), clf
        plot(ens_DN-en0,z_DN,'.r-', ens_UP_adj-en0,z_UP,'.g:')
        hold on, zoom on, grid on, axis ij
        title(['UP(gn)--DN(rd) ens\_no = ' int2str(EnsDNtoUP) ', begin yd=' num2str(yd_DN(1))])
        xlabel(['Ensemble Number (DN) -- ' int2str(en0)])
        ylabel(['ADCPs depth / m'])
    end

    pause
    dx = int2str(EnsDNtoUP);
    xx = inputdlg({'y=done,NaN=async'},...
        '(UP-DN) ens\_no offset',1,{dx});
    if ~isempty(xx) && strncmpi(xx{1},'y',1)
        EnCont = 0;
    else
        if ~isempty(xx)
            if strcmp(xx{1},'NaN')
                xx=NaN;
            else
                xx = str2num(xx{1});
                if isempty(xx)
                    disp('invalid entry, try again after this re-plot')
                end
            end
        end
        if isempty(xx)
            xx = EnsDNtoUP; % either accepted default, or invalid entry
        end
        if ~isnan(xx)
            xx = round(xx);
            if ens_UP(1)-xx>=ens_DN(end) || ens_UP(end)-xx<=ens_DN(1)
                disp('Too large, try again after replotting ...')
            else
                EnsDNtoUP = xx;
            end
            EnCont = 1;
            ens_DN = ens_DN_sav; % in case original were not monotonic (to re-do)
        else
            disp('Continue with asynchronous ADDN and ADUP:')
            disp('Must find CTD time offsets for each separately.')
            EnsDNtoUP = NaN;
            DelDNtoUP = NaN;
            EnCont = 0;
        end
    end % of checking for new ens_offset, or y=done
end % of plotting to check UP-DN synch (by ens_no)
figure(Fens), zoom off
if Fzad
    figure(Fzad), zoom off
end
figure(Fxtra); closereq;

if ~isnan(EnsDNtoUP)
    % replace with adjusted ens_no, offset to sync with ADDN
    ens_UP = ens_UP_adj;
    % initialize to coincide with ADDN time (via sync'd ensemble numbers),
    % offset for ping sequence timing
    ig = find(~isnan(ens_DN+yd_DN));
    yd_UP = interp1(ens_DN(ig), yd_DN(ig), ens_UP) + DelDNtoUP;
end

%% Okay, now get CTD fall rate by 1-sec diff of 24-Hz pressures
SP=get_SWIMS_SciData(yd_b-ex_yday, yd_e+ex_yday, ...
    ['C:\swims\' cruise '\indexes\CTD_' cruise '_matfiles.mat'], ...
    ['C:\swims\' cruise '\data_mat\CTD'], 1);
clear SW
SW.Pr = medfilt1(SP.Pr,25);
SW.yday = SP.yday_adj;
SW.Pitch = SP.Pitch;
SW.Roll = SP.Roll;
SW.Temp = SP.T1;
SW.Sal = SP.S1;
clear SP
SW.Svel = medfilt1(sw_svel(SW.Sal,SW.Temp,SW.Pr),25);
ptsC = 24; % number of 24-Hz points for center difn fall rate (here = 1 s)
%SW.Wpr = 100*(SW.Pr(ptsC+1:end)-SW.Pr(1:end-ptsC)) / (ptsC/24); % Mpa
SW.Wpr = 1*(SW.Pr(ptsC+1:end)-SW.Pr(1:end-ptsC)) / (ptsC/24); % dbar
SW.Wdy = (SW.yday(ptsC+1:end)+SW.yday(1:end-ptsC)) / 2;
ix = find(SW.Wdy<yd_b | SW.Wdy>yd_e);
SW.Wpr(ix) = NaN; SW.Wdy(ix) = NaN;

%% Offset in seconds to add to ADDN (and ADUP) to match SWIMS (PC time)
DelDNtoSW = 0; % seconds (CTD-ADDN)
DelUPtoSW = 0; % seconds, will match DelDNtoSW unless async'd
YD_0 = min(SW.Wdy);

SyCont = 1;
while SyCont
    figure(Fsyn), clf
    sPR = (SW.Wdy-YD_0)*86400;
    sPz = (SW.yday-YD_0)*86400; % for pr vs time
    sDN = (yd_DN-YD_0)*86400 + DelDNtoSW;
    sBT = (yd_BT-YD_0)*86400 + DelDNtoSW;
    sUP = (yd_UP-YD_0)*86400 + DelUPtoSW;
    if isnan(EnsDNtoUP)
        plot(sUP,W_UP,'c-'), hold on
        ww = [W_DN W_BT];
        ss = [sDN sBT];
        ix = find(~isnan(ww+ss));
        ww = ww(ix); ss = ss(ix);
        [sw,iw] = sort(ss);
        plot(sw, ww(iw), 'm-')
    else
        ww = [W_DN W_BT W_UP];
        ss = [sDN sBT sUP];
        ix = find(~isnan(ww+ss));
        ww = ww(ix); ss = ss(ix);
        [sw,iw] = sort(ss);
        plot(sw, ww(iw), 'r-'), hold on
    end
        
    plot(sDN,W_DN,'r.', sBT,W_BT,'k.', sUP,W_UP,'g.'),hold on
    plot(sPR,SW.Wpr,'b-'),grid on,zoom on
    title(['CTD ydays starts ' num2str(YD_0) ', CTD-ADDN = ' num2str(DelDNtoSW) ' s'])
    xlabel('Elapsed seconds (DN=rd,BT=bk,UP=gn, CTD=blue)')
    ylabel('W / m/s')
    
    if Fzad % can plot ADCP-vs-CTD pressures also
        figure(Fzad), clf
        if isnan(EnsDNtoUP)
            plot(sUP,z_UP,'c-', sDN,z_DN,'m-'), hold on, axis ij
        else
            ww = [z_DN z_UP];
            ss = [sDN sUP];
            ix = find(~isnan(ww+ss));
            ww = ww(ix); ss = ss(ix);
            [sw,iw] = sort(ss);
            plot(sw, ww(iw), 'r-'), hold on, axis ij
        end
        plot(sDN,z_DN,'.r', sUP,z_UP,'.g')
        %plot(sPz, SW.Pr*100,'b-') % Mpa
        plot(sPz, SW.Pr,'b-') % dbar
        hold on, zoom on,grid on, axis ij
        title(['CTD ydays starts ' num2str(YD_0) ', CTD-ADDN = ' num2str(DelDNtoSW) ' s'])
        xlabel('Elapsed seconds (DN=rd,BT=bk,UP=gn, CTD=blue)')
        ylabel(['depth(Pr) / m'])
    end

    xd = num2str(DelDNtoSW);
    xu = num2str(DelUPtoSW);
    pause
    if ~isnan(EnsDNtoUP)
        pmts = {'CTD-DN  [y=done]'};
        dfts = {xd};
    else % a-sync mode
        pmts = {'CTD-DN  [y=done]', 'CTD-UP'};
        dfts = {xd, xu};
    end
    xx = inputdlg(pmts, 'CTD sync offsets / sec', 1, dfts);
    if isempty(xx) || strncmpi(xx{1},'y',1)
        SyCont = 0;
    else
        xxd = str2num(xx{1});
        if isempty(xxd)
            disp('invalid entry for CTD-DN, replot ... ')
        elseif abs(xxd)>120
            disp('CTD-DN cannot be more than 120 secs, replot ... ')
        else
            DelDNtoSW = xxd;
        end
        if isnan(EnsDNtoUP) % a-sync mode
            xxu = str2num(xx{2});
            if isempty(xxu)
                disp('invalid entry for CTD-UP, replot ... ')
            elseif abs(xxu)>240
                disp('CTD-UP cannot be more than 240 secs, replot ... ')
            else
                DelUPtoSW = xxu;
            end
        else
            DelUPtoSW = DelDNtoSW; % UP sync'd to DN already
        end
        SyCont = 1;
    end % of checking for new time offset(s), or y=done

end % of plotting to check CTD-DN,CTD-UP synch (by time)
figure(Fens),figure(Fsyn),figure(Fxtra),close([Fens Fsyn Fxtra])
if Fzad
    figure(Fzad),close(Fzad)
end