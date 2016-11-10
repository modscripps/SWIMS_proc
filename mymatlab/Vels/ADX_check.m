% ADX_check.m - Apply error/quality checks to SWIMS ADCP, for BT and
% profile data already in instrument coordinates and in grid matrices

%%% Set bad/error parameters
%% Bottom Track measurement thresholds/limits
Thr_werrBT = 0.2; % error velocity (RDI = BE 1000 mm/s)
Thr_wBT = 4; % vertical speed (cycling)
Thr_uvBT = 5; % horizontal speed 
%% Profile data thresholds/limits:  measured velocities 
%%   should be limited by RDI ambiquity setting = WV 330 to 400 cm/s
Thr_werr = 0.60 ; % error velocity (RDI = WE 2000 mm/s)
Thr_w = 3; % vertical speed (cycling)
Thr_uv = 5; % horizontal speed -- measured velocities
% Depth extents:
Thr_Sfc = 11; % near surface contaminated
Thr_BotPct = 0.92; % sidelobes contaminate last 8% of range-to-bottom
Thr_BotExd = 4; % number of extra meters above Thr_BotPct to exclude
Thr_BinMin = 3; % minimum number of non-NaN bins for various error checks

%% Consistency check:  check for deviations in each w-profile w(:), to
%% catch regions of interference(?) seen in echo contours.
%% Compare w(:) to SWIMS pres.rate (W_pr) for reasonableness, find outliers
Con_wTim = 3 / 86400; % Compare W_pr's within ADCP.yday +/- Con_wTim
Con_wOff = 0.15; % Allow w(:) of W_pr +/- Con_wOff (m/s)
Con_wStd0 = 0.20; % Minimum allowable deviation from median(w(:))
Con_wStds = 2; % number of std dev's allowed from median(w(:))
Con_eRng = 48; % for home02: hi-echos in 54-64-m range, look at>=48m from ADDN
Con_eAlt = 165; % interference occurs once altitude is > 165-m
Con_eAmp = 140; % echo intensities above this are likely bad (if further than Con_eRng)

%% Status error bits, to flag problems before changing to NaNs
Sta_werr = 1;
Sta_w = 2;
Sta_uv = 4;
Sta_Lim = 8;
Sta_wCon = 16;
Sta_eAmp = 32;

% EcAR = [175,25]; % echo>A,dist>R = bad
% difPtRl_Thresh = 4; % badflag pings where pitch/roll changed by more, ping-to-ping
% 
% ix = find( abs(diff(ADX.pitch_DN)) > difPtRl_Thresh | ...
%     abs(diff(ADX.roll_DN)) > difPtRl_Thresh | ...
%     abs(diff(ADX.pitch_UP)) > difPtRl_Thresh | ...
%     abs(diff(ADX.roll_UP)) > difPtRl_Thresh );

%% First, Screen Bottom Track Data
ix = find( abs(ADX.werrBT) > Thr_werrBT | abs(ADX.wBT) > Thr_wBT | ...
    abs(ADX.uBT) > Thr_uvBT | abs(ADX.vBT) > Thr_uvBT );
ADX.werrBT(ix) = NaN; ADX.wBT(ix) = NaN;
ADX.uBT(ix) = NaN; ADX.vBT(ix) = NaN;
%% Screen profile data, but just flag via status bits at first:
ADX.status = 0 * (ADX.u+ADX.v);
% error vel
ix = find( abs(ADX.werr(:)) > Thr_werr );
ADX.status(ix) = ADX.status(ix) + Sta_werr;
% extreme w's
ix = find( abs(ADX.w(:)) > Thr_w );
ADX.status(ix) = ADX.status(ix) + Sta_w;
% extreme u's or v's
ix = find( abs(ADX.u(:)) > Thr_uv | abs(ADX.v(:)) > Thr_uv );
ADX.status(ix) = ADX.status(ix) + Sta_uv;
% Consistency check via vertical velocity: this can reveal poor measurements
% associated with grossly elevated echo intensities;
% Also, flag contaminated bins near surface and near bottom
for icol = 1:length(ADX.yday_DN)
%     zGup = ADX.pr_UP(icol)*100 - 2; % first UP bin above this
%     zGdn = ADX.pr_DN(icol)*100 + 2; % first DN bin below this
    zGup = ADX.pr_UP(icol) - 2; % first UP bin above this
    zGdn = ADX.pr_DN(icol) + 2; % first DN bin below this

    % contaminated bins, sfc/bottom
%     rbot = ADX.bottomBT(icol) - (ADX.pr_BT(icol)*100); % range
%     zmx = ADX.pr_DN(icol)*100 + (rbot*Thr_BotPct) - Thr_BotExd; % contaminated below
    rbot = ADX.bottomBT(icol) - (ADX.pr_BT(icol)); % range
    zmx = ADX.pr_DN(icol) + (rbot*Thr_BotPct) - Thr_BotExd; % contaminated below

    ix = find( (ADX.depth<zGup & ADX.depth<Thr_Sfc) | ...
        (ADX.depth>zGdn & ADX.depth>=zmx) );
    ADX.status(ix,icol) = ADX.status(ix,icol) + Sta_Lim;
    % Check UP/DN bins separately for consistent w's
    for ud=1:2
        % skip 'empty' profiles (already all NaN's)
        if ud==1 & isnan(zGup)
            continue
        elseif ud==2 & isnan(zGdn)
            continue
        end
        switch ud
            case 1 % UP, check non-NaN,non-flagged bins
                izs = find( ADX.depth<zGup & ADX.status(:,icol)==0 );
                izs = flipud(izs); % for UP, to start closest to SWIMS
                iecs = find( ADX.depth<zGup-Con_eRng & ~isnan(ADX.ec1(:,icol)+...
                    ADX.ec2(:,icol)+ADX.ec3(:,icol)+ADX.ec4(:,icol)) );
                iecs = flipud(iecs); % for UP, to start closest to SWIMS
                iecs = []; % No check for UP, for now
                % find SWIMS w=dpr/dt close to profile time
                iw = find( ADX.yday_UP >= ADX.yday_UP(icol)-Con_wTim & ...
                    ADX.yday_UP <= ADX.yday_UP(icol)+Con_wTim );
                wmin = min(ADX.w_UP(iw)) - Con_wOff;
                wmax = max(ADX.w_UP(iw)) + Con_wOff;
            case 2 % DN, check non-NaN,non-flagged bins
                izs = find( ADX.depth>zGdn & ADX.status(:,icol)==0 );
                iecs = [];
                % potential bins for interference, if altitude is enough
                %  if ADX.bottomBT(icol)-ADX.pr_BT(icol)*100>Con_eAlt
                if ADX.bottomBT(icol)-ADX.pr_BT(icol)>Con_eAlt
                    iecs = find( ADX.depth>zGdn+Con_eRng & ~isnan(ADX.ec1(:,icol)+ ...
                        ADX.ec2(:,icol)+ADX.ec3(:,icol)+ADX.ec4(:,icol)) );
                end
                % find SWIMS w=dpr/dt close to profile time
                iw = find( ADX.yday_DN >= ADX.yday_DN(icol)-Con_wTim & ...
                    ADX.yday_DN <= ADX.yday_DN(icol)+Con_wTim );
                wmin = min(ADX.w_DN(iw)) - Con_wOff;
                wmax = max(ADX.w_DN(iw)) + Con_wOff;
        end
        % check w(:) profile for nominal consistency (relative to SWIMS rate)
        iw = []; iwL = 0; wmed=NaN;
        iX = ones(size(izs));  % 1=flag, 0=pass
        if ~isempty(izs)
            iw = find( ADX.w(izs,icol)>wmin & ADX.w(izs,icol)<wmax );
        end
        iwL = length(iw);
        % Get median value of reasonable w(:)'s
        if iwL == 1
            wmed = ADX.w(izs(iw(1)),icol);
        elseif iwL > 1 
            wmed = median( ADX.w(izs(iw), icol) ); % for all ok w(:)
            wmed0 = median( ADX.w(izs(iw(1:min(4,iwL))), icol) ); % nearest to SWIMS
            if abs(wmed-wmed0) > Con_wStd0
                wmed = wmed0; % when inconsistent, use nearest values
            end
        end
        if ~isnan(wmed) % Now, check for outliers within entire w-profile
            % eliminate gross outliers from stats
            iw = find( abs(ADX.w(izs,icol)-wmed) < (Con_wStds+.5)*Con_wStd0 );
            iwL = length(iw);
            % check those remaining
            dw = ADX.w(izs(iw),icol) - wmed; % deviations from median
            sdw = sqrt( sum(dw.^2) / iwL ); % std dev, sort of
            sdw = max(Con_wStds * sdw, Con_wStd0); % acceptable deviations
            iok = find(abs(dw) <= sdw); % okay w(:)s
            if length(iok)<iwL % exclude outliers to recompute tolerance
                sdw0 = sqrt( sum(dw(iok).^2) / length(iok) ); % new std dev
                sdw0 = max(Con_wStds * sdw0, Con_wStd0);
                sdw = min(sdw,sdw0); % new acceptable deviation
            end
            iok = find(abs(dw) <= sdw);
            iX( iw(iok) ) = 0; % These passed the consistency check
        end
        ix = find( iX>0 );
        % Flag outliers (UP or DN profile)
        if ~isempty(ix)
            ADX.status(izs(ix),icol) = ADX.status(izs(ix),icol) + Sta_wCon;
        end
        %
        %% Next, check for hi-echo patches below ADDN
%         if ~isempty(iecs)
%             % find bins with at least one beam of high echo intensity
%             ie = find( max( [ADX.ec1(iecs,icol) ADX.ec2(iecs,icol) ...
%                     ADX.ec3(iecs,icol) ADX.ec4(iecs,icol)], [], 2 ) > Con_eAmp);
%             if ~isempty(ie)
%                 % find 'Patches' = two or more consecutive bins
%                 ieP = find(diff(iecs(ie))==1); % 'Patches' = two or more consecutive
%                 if ~isempty(ieP)
%                     % flag bins adjacent to patches also
%                     ie = [ie; ie(ieP)-1; ie(ieP+1)+1]; 
%                 end
%                 ix = unique(ie);
%                 ie = find(ix<1 | ix>length(iecs));
%                 ix(ie) = [];  % if 'adjacent' to first or last bin
%                 ADX.status(iecs(ix),icol) = ADX.status(iecs(ix),icol) + Sta_eAmp;
%             end
%         end
    end % of ud=1:2 loop; UP then DN profile
end % of icol-loop; profiles
                
