% ADX_Prof_avg_ToRef.m - with single ping velocities in earth coords, average
% into profiles, using supplied profile REF.{u,v,depth_range, Z_ok(range)} as a seed

% initialize
ADX.uAbs = NaN*ADX.u; ADX.vAbs = NaN*ADX.v; % absolute velocities (using REF)
ADX.uRel = NaN*ADX.u; ADX.vRel = NaN*ADX.v; % relative results, with unknown offset
% Also, compute mean on bin-to-bin diffs, to integrate shear as 2nd refn profile
ADX.du = diff(ADX.u);
ADX.dv = diff(ADX.v);

% Determine whether SWIMS profile is upward or downward
igd = find(~isnan(ADX.pr_DN));
igu = find(~isnan(ADX.pr_UP));
if length(igd)<4
    if length(igu)<4
        disp('Not enough pings to average')
        return
    else
        xx = ADX.pr_UP(igu(end))-ADX.pr_UP(igu(1));
    end
else
    xx = ADX.pr_DN(igd(end))-ADX.pr_DN(igd(1));
end
if xx > 0
    UpDown = 2; % Down profile
else
    UpDown = 1; % Up profile
end

% Record bin numbers just above ADUP and below ADDN for each ensemble
fBup = 1 * ones(size(ADX.yday_UP)); % default to shallowest bin
fBdn = length(ADX.depth) * ones(size(ADX.yday_DN)); % default to deepest bin
for ip = 1:length(ADX.yday_BT)
    % zGup = ADX.pr_UP(ip)*100 - 2; % first UP bin is above this
    zGup = ADX.pr_UP(ip) - 2;
    nb = find(ADX.depth < zGup); 
    if ~isempty(nb)
        fBup(ip) = nb(end); % first possible bin above ADUP
    end
    % zGdn = ADX.pr_DN(ip)*100 + 2; % first DN bin is below this
    zGdn = ADX.pr_DN(ip) + 2;
    nb = find(ADX.depth > zGdn); 
    if ~isempty(nb)
        fBdn(ip) = nb(1); % first possible bin below ADDN
    end
end % of determining starting bins

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Map suitable interval of reference profiles onto output grid
uAe = 0 * ADX.depth;
vAe = uAe;
ig = find( ADX.depth >= REF.Z_ok(1) & ADX.depth <= REF.Z_ok(2) );
if length(ig)<2
    AVG=[];
    return
end
ik = find(~isnan(REF.u+REF.v)); % just in case, interp over NaNs
uAe(ig) = interp1(REF.depth(ik), REF.u(ik), ADX.depth(ig));
vAe(ig) = interp1(REF.depth(ik), REF.v(ik), ADX.depth(ig));

nMin = 8; % Minimum number of values accumulated in bin to use in average
uShe = 0 * ADX.depth; vShe = uShe; % shear estimates
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
AVG.U_rel2abs = NaN; % REF-based offset, if applied
AVG.V_rel2abs = NaN; % REF-based offset, if applied
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

%% Compute mean profile of shears from all ensembles (for intermediate estimate)
for ib=1:length(ADX.depth)
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
end % of looping thru bins, computing diffs,shears
% integrate diffs(uv) to form mean velocity profiles with unknown offsets
uRe = cumsum(uShe); vRe = cumsum(vShe);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Use mid-depth interval of these profiles as initial reference;
% Compute adjustments for all profiles by matching weighted averages

%% First, based on SWIMS ADCP shear-integrated profiles:
ig = find(bCsh>9); iC = round( (ig(1)+ig(end)) / 2 );
ig = min( round( (ig(end)-ig(1))/4 ), 20) ; % middle half, but +/-20 bins at most
iSR = [iC-ig:iC+ig]';
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

%% Second, based on supplied REFerence profiles (ship based, often, and presumed
%%  to be absolute):
iRF = find( ADX.depth >= REF.Z_ok(1) & ADX.depth <= REF.Z_ok(2) );
iC = round( (iRF(1)+iRF(end)) / 2 ); % center bin of refn interval
% Set counts to use only this interval initially
bCt = 0*bCt; bCt(iRF) = nMin*2; % initialize bins with equal weights
% find various groups of ensembles, relative to initial reference interval
iSend = length(ADX.yday_UP);
% iSabv = find(ADX.pr_UP*100<min(ADX.depth(iRF)) ); % above interval
iSabv = find(ADX.pr_UP<min(ADX.depth(iRF)) );
if isempty(iSabv)
    iSabv = [iSend 1]; iSabv = iSabv(UpDown); % shallowest
end
% iSbel = find(ADX.pr_DN*100>max(ADX.depth(iRF)) ); % below interval
iSbel = find(ADX.pr_DN>max(ADX.depth(iRF)) );
if isempty(iSbel)
    iSbel = [1 iSend]; iSbel = iSbel(UpDown); % deepest
end
% [x, iSmid] = min( abs(ADX.pr_DN*100-ADX.depth(iC)) );
[x, iSmid] = min( abs(ADX.pr_DN-ADX.depth(iC)) );
% Compute in order: DN in upper half of interval, UP in lower half,
%   DN above interval upward, UP above mid-int upward,
%   UP below interval downward, DN below mid-int downward (1=UP,2=DN)
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
dup = [2, 1, 2, 1, 1, 2]; % directions
%% Accumulate profile #s, directions (after shear-seed based ones)
for i=1:6
    if ~isempty(ii{i})
        DoProfs = [DoProfs, [ii{i}; dup(i)*ones( size(ii{i}) ) ] ];
    end
end

%keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Determine offsets for UP and DN profiles by matching weighted
%% averages to estimated mean profiles:
%% UP vs REF, DN vs REF, UP vs integ(Sh), DN vs integ(Sh)
uOff = NaN * [ADX.yday_DN; ADX.yday_DN; ADX.yday_DN; ADX.yday_DN]; 
vOff = uOff;
WTnRF = 1; % as accumulating, add this much weight for each SWIMS ensemble
WTvSH = 8/5; % catches up with integ(Sh) after 5 profiles
% As offsets are estimated, accumulate adjusted profiles into average profiles
for ip=1:length(DoProfs(1,:))
    iprof = DoProfs(1,ip); % column no.
    tyProf = DoProfs(2,ip); % type
    switch tyProf
        case 1 % ADUP, vs REF
            iG = find( ADX.depth<=ADX.depth(fBup(iprof)) & ...
                ~isnan(ADX.u(:,iprof)+ADX.v(:,iprof)) & bCt>=nMin ); % bins for matching
            iGb = find( ADX.depth<=ADX.depth(fBup(iprof)) & ...
                ~isnan(ADX.u(:,iprof)+ADX.v(:,iprof)) ); % bins to adjust
            namU = 'uAe'; namV = 'vAe'; namCt = 'bCt'; WtProf = WTnRF; namuv='Abs';
        case 2 % ADDN, vs REF
            iG = find( ADX.depth>=ADX.depth(fBdn(iprof)) & ...
                ~isnan(ADX.u(:,iprof)+ADX.v(:,iprof)) & bCt>=nMin ); % bins for matching
            iGb = find( ADX.depth>=ADX.depth(fBdn(iprof)) & ...
                ~isnan(ADX.u(:,iprof)+ADX.v(:,iprof)) ); % bins to adjust
            namU = 'uAe'; namV = 'vAe'; namCt = 'bCt'; WtProf = WTnRF; namuv='Abs';
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
%%%%%%%%%%%%%%%%%%%% Final Averaging %%%%%%%%%%%%%%%%%%%%%%%%%%
%% Now, Re-Average adjusted profiles, rejecting outliers
%% and saving bin stats
for ib=1:length(AVG.depth)
    for ity = 1:2 % absolute, relative profiles
        if ity==1
            uu = ADX.uAbs(ib,:); vv = ADX.vAbs(ib,:); % offset via REF
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Offset relative profiles based on interval for REF'd profiles
ig = find( AVG.depth >= REF.Z_ok(1) & AVG.depth <= REF.Z_ok(2) & ...
    ~isnan(AVG.U_abs+AVG.V_abs+AVG.U_rel+AVG.V_rel) );
% Compute, apply, and save offsets to make relative profiles absolute
if length(ig) > nMin/2 % Only if good bins are available
    uRef = mean(AVG.U_abs(ig));
    uPro = mean(AVG.U_rel(ig));
    vRef = mean(AVG.V_abs(ig));
    vPro = mean(AVG.V_rel(ig));
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
    

                    
            