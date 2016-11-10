% ADX_Prof_avg.m - with single ping velocities in earth coords, avg into profiles

clear ADDN ADUP xx yy
ADX.w12=[]; ADX.w43=[]; ADX.w=[]; ADX.werr=[]; % clear some room
% absolute velocities 
ig = find(~isnan(ADX.yday_BT));
uBTup = interp1(ADX.yday_BT(ig), ADX.uBT(ig), ADX.yday_UP );
vBTup = interp1(ADX.yday_BT(ig), ADX.vBT(ig), ADX.yday_UP );
uBTdn = interp1(ADX.yday_BT(ig), ADX.uBT(ig), ADX.yday_DN );
vBTdn = interp1(ADX.yday_BT(ig), ADX.vBT(ig), ADX.yday_DN );
ix = find(isnan(ADX.bottomBT));
uBTup(ix)=NaN; vBTup(ix)=NaN;
uBTdn(ix)=NaN; vBTdn(ix)=NaN;
% remove those flagged for errors, unreliable, etc.
ix = find( ~ (ADX.status(:)==0) );
ADX.u(ix) = NaN; ADX.v(ix) = NaN;
% initialize
ADX.uAbs = NaN*ADX.u; ADX.vAbs = NaN*ADX.v; % absolute velocities
ADX.uRel = NaN*ADX.u; ADX.vRel = NaN*ADX.v; % relative results, with unknown offset
% Also, compute mean on bin-to-bin diffs, to integrate where BT fails
ADX.du = diff(ADX.u);
ADX.dv = diff(ADX.v);

if ADX.pr_DN(end)>ADX.pr_DN(1)
    UpDown = 2; % Down profile
else
    UpDown = 1; % Up profile
end

%% Compute absolute velocity for BT'd profiles.
% record bin numbers just above ADUP and below ADDN for each ensemble
fBup = 1 * ones(size(ADX.yday_UP)); % default to shallowest bin
fBdn = length(ADX.depth) * ones(size(ADX.yday_DN)); % default to deepest bin
for ip = 1:length(ADX.yday_BT)
    % remove BT vel from ADUP
    % zGup = ADX.pr_UP(ip)*100 - 2; % first UP bin is above this
    zGup = ADX.pr_UP(ip) - 2;
    nb = find(ADX.depth < zGup); 
    if ~isempty(nb)
        fBup(ip) = nb(end); % first possible bin above ADUP
        % apply offset to bins above ADUP
        ADX.uAbs(1:fBup(ip),ip) = ADX.u(1:fBup(ip),ip) - uBTup(ip);
        ADX.vAbs(1:fBup(ip),ip) = ADX.v(1:fBup(ip),ip) - vBTup(ip);
    end
    % remove BT vel from ADDN
    % zGdn = ADX.pr_DN(ip)*100 + 2; % first DN bin is below this
    zGdn = ADX.pr_DN(ip) + 2;
    nb = find(ADX.depth > zGdn); 
    if ~isempty(nb)
        fBdn(ip) = nb(1); % first possible bin below ADDN
        % apply offset to bins below ADDN
        ADX.uAbs(fBdn(ip):end,ip) = ADX.u(fBdn(ip):end,ip) - uBTdn(ip);
        ADX.vAbs(fBdn(ip):end,ip) = ADX.v(fBdn(ip):end,ip) - vBTdn(ip);
    end
end % of removing BT velocities

nMin = 8; % Minimum number of values accumulated in bin to use in average
uAe = 0 * ADX.depth;
vAe = uAe;
uShe = uAe; vShe = uAe; % shear estimates
bCt = 0 * ADX.depth;
bCsh = bCt;
%% initialize result structure AVG
clear AVG
AVG.depth = ADX.depth;
AVG.yd_min = NaN*AVG.depth;
AVG.yd_max = AVG.yd_min;
AVG.U_abs = NaN*AVG.depth;
AVG.V_abs = AVG.U_abs;
AVG.U_rel = AVG.U_abs;
AVG.V_rel = AVG.U_abs;
AVG.U_abs_std = AVG.U_abs;
AVG.V_abs_std = AVG.U_abs;
AVG.U_rel_std = AVG.U_abs;
AVG.V_rel_std = AVG.U_abs;
AVG.count_abs = AVG.U_abs;
AVG.count_rel = AVG.U_abs;
AVG.outlier_abs = AVG.U_abs;
AVG.outlier_rel = AVG.U_abs;
AVG.U_rel2abs = NaN; % BT-based offset, if applied
AVG.V_rel2abs = NaN; % BT-based offset, if applied
nSTD = 2.3; % number of std deviations for outliers (velocity)
% For averaged shear, set default=0 for later integration,
% but stderr,count = NaN to show no data were found
AVG.depSH = (AVG.depth(1:end-1)+AVG.depth(2:end)) / 2;
AVG.dUdz = 0 * ones(size(AVG.depSH));
AVG.dVdz = AVG.dUdz;
AVG.dU_std = NaN*AVG.dUdz;
AVG.dV_std = AVG.dU_std;
AVG.count_dUV = AVG.dU_std;
AVG.outlier_dUV = AVG.dU_std;
dZ = median(diff(AVG.depth)); % bin size (m)
nSTD_sh = 2.0; % nSTD (for shear)
%% Compute estimate of average absolute velocity using BT'd values;
% Also, compute mean profile of shears from all ensembles 
for ib=1:length(ADX.depth)
    ig = find( ~isnan(ADX.uAbs(ib,:)+ADX.vAbs(ib,:)) );
    bCt(ib) = length(ig);
    if ~isempty(ig)
        uAe(ib) = mean(ADX.uAbs(ib,ig));
        vAe(ib) = mean(ADX.vAbs(ib,ig));
    end
    % take mean of shears, but exclude 20% highest and lowest
    % (first bin stays=0 for cumsum). Re-avg for final avgd shear profiles,
    % but use integral of intermediate results as 'seed' for estimating
    % ping-offsets for relative velocity calculation.
    if ib<=length(ADX.du(:,1))
        ig = find( ~isnan(ADX.du(ib,:)+ADX.dv(ib,:)) );
        bCsh(ib+1) = length(ig);
        if bCsh(ib+1)>9
            uShe(ib+1) = trimmean(ADX.du(ib,ig), 40);
            vShe(ib+1) = trimmean(ADX.dv(ib,ig), 40);
        elseif bCsh(ib+1)>0
            uShe(ib+1) = mean(ADX.du(ib,ig));
            vShe(ib+1) = mean(ADX.dv(ib,ig));
        end
        %% compute std dev from calc'd mean, exclude outliers, re-avg for
        %% final shear profiles ( in units = 1/s )
        if bCsh(ib+1) >= max(nMin, 3)
            du = ADX.du(ib,ig) - uShe(ib+1);
            dv = ADX.dv(ib,ig) - vShe(ib+1);
            stu = sqrt( sum(du.^2) / (bCsh(ib+1)-1) ); % like std, type=1
            stv = sqrt( sum(dv.^2) / (bCsh(ib+1)-1) );
            % determine outliers (exclude if in either u,v component)
            iex = find( abs(du) > nSTD_sh*stu | abs(dv) > nSTD_sh*stv );
            cto = length(iex);
            if cto>0 % outliers were found, exclude for re-calc
                ig(iex) = [];
            end
            % recalc using all remaining (result = shear, not just diff)
            AVG.dUdz(ib) = mean(ADX.du(ib,ig) / dZ);
            AVG.dVdz(ib) = mean(ADX.dv(ib,ig) / dZ);
            if length(ig)>2
                AVG.dU_std(ib) = std(ADX.du(ib,ig)/dZ, 1);
                AVG.dV_std(ib) = std(ADX.dv(ib,ig)/dZ, 1); 
            end
            AVG.count_dUV(ib) = bCsh(ib+1);
            AVG.outlier_dUV(ib) = cto;
        elseif bCsh(ib+1) > 0
            % record that some were found, but none were used
            AVG.count_dUV(ib) = bCsh(ib+1);
            AVG.outlier_dUV(ib) = bCsh(ib+1);
        end % of portion where final avgd shears are computed
    end % of portion where avg'd shears or diffs(vel) are computed    
end % of looping thru bins, computing abs vels or diffs,shears
% integrate diffs(uv) to form mean velocity profiles with unknown offsets
uRe = cumsum(uShe); vRe = cumsum(vShe);

%% Use mid-depth interval of these profiles as initial reference;
% Compute adjustments for all profiles by matching weighted averages
ig = find(bCsh>9); 
if ~isempty(ig)
    iC = round( (ig(1)+ig(end)) / 2 );
    ig = min( round( (ig(end)-ig(1))/4 ), 20) ; % middle half, but +/-20 bins at most
    iSR = [iC-ig:iC+ig]';
else disp('ig is empty')
    iSR = 0;
end

% Reset counts,vels to use only this interval initially
uSRe = 0*ones(size(uRe)); uSRe(iSR) = uRe(iSR);
vSRe = 0*ones(size(vRe)); vSRe(iSR) = vRe(iSR);
bCsh = 0*bCsh; bCsh(iSR) = nMin; % initialize bins with equal weights
% find various groups of ensembles, relative to initial reference interval
iSend = length(ADX.yday_UP);
% iSabv = find(ADX.pr_UP*100<min(ADX.depth(iSR)) ); % above interval
iSabv = find(ADX.pr_UP<min(ADX.depth(iSR)) );
if isempty(iSabv)
    iSabv = [iSend 1]; iSabv = iSabv(UpDown); % shallowest
end
% iSbel = find(ADX.pr_DN*100>max(ADX.depth(iSR)) ); % below interval
iSbel = find(ADX.pr_DN>max(ADX.depth(iSR)) );
if isempty(iSbel)
    iSbel = [1 iSend]; iSbel = iSbel(UpDown); % deepest
end
% [x, iSmid] = min( abs(ADX.pr_DN*100-ADX.depth(iC)) );
[x, iSmid] = min( abs(ADX.pr_DN-ADX.depth(iC)) );
    
% Compute in order: DN in upper half of interval, UP in lower half,
%   DN above interval upward, UP above mid-int upward,
%   UP below interval downward, DN below mid-int downward (3=UP,4=DN)
clear ii
if UpDown==1 % Upward profile
    iSabv = iSabv(1); iSbel = iSbel(end); % bins abutting interval
    ii{1} = iSabv-1:-1:iSmid; ii{2} = iSbel+1:iSmid;
    ii{3} = iSabv:iSend; ii{4} = iSmid+1:iSend;
    ii{5} = iSbel-1:-1:1; ii{6} = iSmid-1:-1:1;
else % Downward profile
    iSabv = iSabv(end); iSbel = iSbel(1); % bins abutting interval
    ii{1} = iSabv+1:iSmid; ii{2} = iSbel-1:-1:iSmid;
    ii{3} = iSabv:-1:1; ii{4} = iSmid-1:-1:1;
    ii{5} = iSbel:iSend; ii{6} = iSmid+1:iSend;
end  
dup = [4, 3, 4, 3, 3, 4]; % directions
DoProfs = [];
%% Accumulate profile #s, directions
for i=1:6
    if ~isempty(ii{i})
        DoProfs = [DoProfs, [ii{i}; dup(i)*ones( size(ii{i}) ) ] ];
    end
end

%% Adjust non-BT'd profiles by matching weighted average of good bins
iBT = find( ~isnan(uBTdn+vBTup) ); % BT'd (okay for both UP and DN)
iNoBT = find( isnan(uBTdn+vBTup) );
uAe_sv = uAe; vAe_sv = vAe; bCt_sv = bCt; % save if needed later
if length(iBT) > nMin & ~isempty(iNoBT) % Make sure there are some to do
    % find non-BT'd within, before, and after good BT interval
    iSmid = find(iNoBT>min(iBT) & iNoBT<max(iBT));
    iSmid = iNoBT(iSmid); % those within interval
    % For marginal BT'ing, re-do even BT'd ones:
    ik = find(bCt>0); ig = find(bCt>=nMin); % total BT bins, useful BT bins
    if length(ik)>nMin & length(ig)<nMin
        ik = find(bCt>0 & bCt<nMin);
        bCt(ik) = nMin; % force minimal useful weight for these
        iSmid = iBT(1):iBT(end); % all ensembles in BT span
    end
    iSend = length(ADX.yday_UP);
    % Compute in order. Within BT interval: DN downward then UP upward;
    %   above interval: DN upward then UP upward,
    %   below interval: UP downward then DN downward (1=UP,2=DN)
    clear ii
    if UpDown==1 % Upward profile
        ii{1} = fliplr(iSmid); ii{2} = iSmid;
        ii{3} = iBT(end)+1:iSend; ii{4} = ii{3};
        ii{5} = iBT(1)-1:-1:1; ii{6} = ii{5};
    else % Downward profile
        ii{1} = iSmid; ii{2} = fliplr(iSmid);
        ii{3} = iBT(1)-1:-1:1; ii{4} = ii{3};
        ii{5} = iBT(end)+1:iSend; ii{6} = ii{5};
    end
    dup = [2, 1, 2, 1, 1, 2]; % directions
    for i=1:6
        if ~isempty(ii{i})
            DoProfs = [DoProfs, [ii{i}; dup(i)*ones( size(ii{i}) ) ] ];
        end
    end
end
       
%keyboard

%% Determine offsets for UP and DN profiles by matching weighted
%% averages to estimated mean profiles:
%% UP vs BT, DN vs BT, UP vs integ(Sh), DN vs integ(Sh)
uOff = [uBTup; uBTdn; NaN*uBTup; NaN*uBTdn]; 
vOff = [vBTup; vBTdn; NaN*vBTup; NaN*vBTdn];
WTnBT = 0.5; % as accumulating, give non-BTs half the weight of (initials) BTs
WTvSH = 8/5; % catches up with integ(Sh) after 5 profiles
% As offsets are estimated, accumulate adjusted profiles into average profiles
for ip=1:length(DoProfs(1,:))
    iprof = DoProfs(1,ip); % column no.
    tyProf = DoProfs(2,ip); % type
    switch tyProf
        case 1 % ADUP, vs BT
            iG = find( ADX.depth<=ADX.depth(fBup(iprof)) & ...
                ~isnan(ADX.u(:,iprof)+ADX.v(:,iprof)) & bCt>=nMin ); % bins for matching
            iGb = find( ADX.depth<=ADX.depth(fBup(iprof)) & ...
                ~isnan(ADX.u(:,iprof)+ADX.v(:,iprof)) ); % bins to adjust
            namU = 'uAe'; namV = 'vAe'; namCt = 'bCt'; WtProf = WTnBT; namuv='Abs';
        case 2 % ADDN, vs BT
            iG = find( ADX.depth>=ADX.depth(fBdn(iprof)) & ...
                ~isnan(ADX.u(:,iprof)+ADX.v(:,iprof)) & bCt>=nMin ); % bins for matching
            iGb = find( ADX.depth>=ADX.depth(fBdn(iprof)) & ...
                ~isnan(ADX.u(:,iprof)+ADX.v(:,iprof)) ); % bins to adjust
            namU = 'uAe'; namV = 'vAe'; namCt = 'bCt'; WtProf = WTnBT; namuv='Abs';
        case 3 % ADUP, vs integ(Sh)
            iG = find( ADX.depth<=ADX.depth(fBup(iprof)) & ...
                ~isnan(ADX.u(:,iprof)+ADX.v(:,iprof)) & bCsh>=nMin ); % bins for matching
            iGb = find( ADX.depth<=ADX.depth(fBup(iprof)) & ...
                ~isnan(ADX.u(:,iprof)+ADX.v(:,iprof)) ); % bins to adjust
            namU = 'uSRe'; namV = 'vSRe'; namCt = 'bCsh'; WtProf = WTvSH; namuv='Rel';
        case 4 % ADDN, vs integ(Sh)
            iG = find( ADX.depth>=ADX.depth(fBdn(iprof)) & ...
                ~isnan(ADX.u(:,iprof)+ADX.v(:,iprof)) & bCsh>=nMin ); % bins for matching
            iGb = find( ADX.depth>=ADX.depth(fBdn(iprof)) & ...
                ~isnan(ADX.u(:,iprof)+ADX.v(:,iprof)) ); % bins to adjust
            namU = 'uSRe'; namV = 'vSRe'; namCt = 'bCsh'; WtProf = WTvSH; namuv='Rel';
    end
    % Copy mean/reference profiles, weights
    eval(['UU = ' namU ';']);
    eval(['VV = ' namV ';']);
    eval(['CC = ' namCt ';']);
    % Compute offsets
    if ~isempty(iG)
        twt = sum(CC(iG)); % total weight
        uRef = sum( UU(iG).*(CC(iG)/twt) );
        uPro = sum( ADX.u(iG,iprof).*(CC(iG)/twt) );
        vRef = sum( VV(iG).*(CC(iG)/twt) );
        vPro = sum( ADX.v(iG,iprof).*(CC(iG)/twt) );
        % infer offsets from weighted averages
        uOff(tyProf, iprof) = uPro - uRef;
        vOff(tyProf, iprof) = vPro - vRef;
    end
    % adjust mean profiles, and total weight per bin; Save adjusted ensembles
    if ~isempty(iGb) & ~isnan( uOff(tyProf, iprof)+vOff(tyProf, iprof) )
        uadj = ADX.u(iGb,iprof)-uOff(tyProf, iprof);
        vadj = ADX.v(iGb,iprof)-vOff(tyProf, iprof);
        UU(iGb) = UU(iGb).*CC(iGb) + WtProf*( uadj );
        VV(iGb) = VV(iGb).*CC(iGb) + WtProf*( vadj );
        CC(iGb) = CC(iGb) + WtProf;
        UU(iGb) = UU(iGb) ./ CC(iGb);
        VV(iGb) = VV(iGb) ./ CC(iGb);
        % copy back into specified vector and array profiles
        eval([namU ' = UU;']);
        eval([namV ' = VV;']);
        eval([namCt ' = CC;']);
        eval(['ADX.u' namuv '(iGb,iprof) = uadj;']);
        eval(['ADX.v' namuv '(iGb,iprof) = vadj;']);
    end
end % of DoProds=[profiles;types] loop for computing offsets 
clear UU VV CC

%% Now, Re-Average adjusted profiles, rejecting outliers
%% and saving bin stats
for ib=1:length(AVG.depth)
    for ity = 1:2 % absolute, relative profiles
        if ity==1 & length(iBT) < 4 % Not enough BT'd profiles, skip
            continue
        end
        if ity==1
            uu = ADX.uAbs(ib,:); vv = ADX.vAbs(ib,:); % offset via BT
        else
            uu = ADX.uRel(ib,:); vv = ADX.vRel(ib,:); % offset via integ(Sh)
        end
        ig = find(~isnan(uu+vv));
        ct = length(ig); cto = 0;
        if isempty(ig) 
            continue; % no data for this bin
        end
        UU = mean(uu(ig)); VV = mean(vv(ig)); US = NaN; VS = NaN;
        % record timespan of pings averaged into this bin (incl/outliers)
        AVG.yd_min(ib) = min([AVG.yd_min(ib), ADX.yday_UP(ig), ADX.yday_DN(ig)]);
        AVG.yd_max(ib) = max([AVG.yd_max(ib), ADX.yday_UP(ig), ADX.yday_DN(ig)]);
        if ct > 2
            US = std(uu(ig),1); VS = std(vv(ig),1);
            % determine outliers (exclude if in either u,v component)
            iex = find( abs(uu(ig)-UU) > nSTD*US | abs(vv(ig)-VV) > nSTD*VS );
            cto = length(iex);
            if cto>0 % outliers were found, recompute mean and std dev
                ig(iex) = [];
                UU = mean(uu(ig)); VV = mean(vv(ig)); US = NaN; VS = NaN;
                if length(ig)>2
                    US = std(uu(ig),1); VS = std(vv(ig),1);
                end
            end
        end
        % Save results for this depth bin
        if ity==1
            AVG.U_abs(ib) = UU;
            AVG.V_abs(ib) = VV;
            AVG.U_abs_std(ib) = US;
            AVG.V_abs_std(ib) = VS;
            AVG.count_abs(ib) = ct;
            AVG.outlier_abs(ib) = cto;
        else
            AVG.U_rel(ib) = UU;
            AVG.V_rel(ib) = VV;
            AVG.U_rel_std(ib) = US;
            AVG.V_rel_std(ib) = VS;
            AVG.count_rel(ib) = ct;
            AVG.outlier_rel(ib) = cto;
        end
    end % of absolute/relative loop
end % of depth bin (ib) loop, for computing final avg profiles

%% Offset relative profiles based on reasonable BT'd profiles
% As a fallback, allow fewer BT'd profiles if not enough bins found
for i_thr=nMin:-1:2
    ig = find(bCt_sv > i_thr & ~isnan(AVG.U_rel+AVG.V_rel) );
    if length(ig) >= nMin
        break; % found enough
    end
end
% if enough bins are available, use those with most corrected data
ik = [];
if length(ig) > 30
    x = sort(bCt_sv(ig));
    ik = find(x > x(1)); if isempty(ik), ik=1; end
    x = max( x(end-15), x(ik(1)) ); % in case of several minimal values
    ik = find(bCt_sv >= x); % use 15-ish bins with most BT'd data
elseif length(ig) > 10
    x = sort(bCt_sv(ig));
    x = max( x(floor(length(ig)/2)), x(1)+1 );
    ik = find( bCt_sv >= x );
end
if length(ik) >= nMin
    ig = ik; % use this subset of bins
end
% Compute, apply, and save offsets to make relative profiles absolute
if length(ig) > nMin/2 % Only if enough BT'd bins are available
    twt = sum(bCt_sv(ig)); % total weight
    uRef = sum( uAe_sv(ig).*(bCt_sv(ig)/twt) );
    uPro = sum( AVG.U_rel(ig).*(bCt_sv(ig)/twt) );
    vRef = sum( vAe_sv(ig).*(bCt_sv(ig)/twt) );
    vPro = sum( AVG.V_rel(ig).*(bCt_sv(ig)/twt) );
    % infer offsets from weighted averages
    AVG.U_rel2abs = uPro - uRef;
    AVG.V_rel2abs = vPro - vRef;
    AVG.U_rel = AVG.U_rel - AVG.U_rel2abs;
    AVG.V_rel = AVG.V_rel - AVG.V_rel2abs;
end
% finally, also save offsets and other info for individual UP/DN profiles
fn = fieldnames(ADX);
AVG.ENSinfo.uOffsets = uOff;
AVG.ENSinfo.vOffsets = vOff;
eLen = size(ADX.yday_DN, 2);
for i=1:length(fn)
    if size(ADX.(fn{i}), 1) == 1 & size(ADX.(fn{i}), 2) == eLen
        AVG.ENSinfo.(fn{i}) = ADX.(fn{i});
    end
end
% clear ADX
AVG.uOffsets = uOff;
AVG.vOffsets = vOff;
    

                    
            