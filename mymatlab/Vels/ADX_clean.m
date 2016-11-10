% ADX_clean.m - Use results of screening/error checks to decide which depth/profile
%  bins of SWIMS ADCP data (in E,N coords) are available for averaging

%% Status error bits, set to flag problems in depth/profile bins
% Existing conditions: (NaN = badflagged by ADCP or gross magnitude checks)
Sta_werr = 1;
Sta_w = 2;
Sta_uv = 4;
Sta_Lim = 8;
Sta_wCon = 16;
Sta_eAmp = 32;
% Additional conditions determined here:
Sta_orphan = 1024; %

%% Parameters for further exclusions:
% Contiguous good bins required, based on range from ADCP:
Rng_sfc = 25; % fewer required if depth<25m
Bins_sfc = [2, 1,inf];
Rng_bot = 25; % fewer required if altitude<25m
Bins_bot = [2, 1,inf];
Rng_Swims = [ [0,40]; [40, inf] ]; % Range (m) of closest bin in chunk
% For each range: [number bins reqd,  min,max index gap before chunk; ... ]
Bins_Reqd = { [3, 1,6; 4, 7,inf]; ...
        [4, 1,1;  5, 2,4;  8, 5,10;  20, 11,inf] };


% Also, flag contaminated bins near surface and near bottom
for icol = 1:length(ADX.yday_DN)
    % find bins with data that passed all tests (ADUP)
%     zGup = ADX.pr_UP(icol)*100 - 2; % first UP bin is above this
    zGup = ADX.pr_UP(icol) - 2; 
    to_sfc = zGup; % range to sfc
    fBup = find(ADX.depth < zGup); 
    if ~isempty(fBup)
        fBup = fBup(end); % first possible bin above ADUP
        iGup = find(ADX.depth<=ADX.depth(fBup) & ADX.status(:,icol)==0 );
    else
        fBup = NaN; iGup = [];
    end
    % find bins with data that passed all tests (ADDN)
    % zGdn = ADX.pr_DN(icol)*100 + 2; % first DN bin is below this
    zGdn = ADX.pr_DN(icol) + 2;
    to_bot = ADX.bottomBT(icol) - zGdn; % range to bottom
    fBdn = find(ADX.depth > zGdn); 
    if ~isempty(fBdn)
        fBdn = fBdn(1); % first possible bin below ADDN
        iGdn = find(ADX.depth>=ADX.depth(fBdn) & ADX.status(:,icol)==0 );
    else
        fBdn = NaN; iGdn = [];
    end
    % Check UP/DN bins separately for contiguous chunks of good data
    for ud=1:2
        % skip 'empty' profiles (already all NaN's or flagged)
        if ud==1 & isempty(iGup)
            continue
        elseif ud==2 & isempty(iGdn)
            continue
        end
        % Use default Gap parameters, unless close to surface or bottom
        Rng_USE = Rng_Swims;
        Bins_USE = Bins_Reqd;
        switch ud
            case 1 % UP, check bins from ADUP upwards
                iG = flipud(iGup); zRef = zGup;
                dgB = diff(iGup); dgB = flipud(dgB); % gaps are > 1
                dgB = [fBup - iGup(end); dgB]; % first one is gap above ADUP
                if to_sfc <= Rng_sfc % close to sfc
                    Rng_USE = [0 Rng_sfc+1];
                    Bins_USE = { Bins_sfc };
                end
            case 2 % DN, check bins from ADDN downwards
                iG = iGdn; zRef = zGdn;
                dgB = diff(iGdn); % gaps are > 1
                dgB = [iGdn(1) - fBdn; dgB]; % first one is gap below ADDN
                if to_bot <= Rng_bot
                    Rng_USE = [0 Rng_bot+1];
                    Bins_USE = { Bins_bot };
                end
        end
        iX = 0 * iG; % default = flag as orphan bins
        LGaps = []; % gap lengths (in bins)
        ICnks = []; % chunks, as indices into iG;
        nc = 0;
        %% loop thru bins, finding gap lengths and chunk indices
        for i = 1:length(iG)
            if dgB(i)~=1 | i==1
                nc = nc+1; % start of new chunk
                ICnks{nc} = i;
                LGaps(nc) = max(dgB(i)-1, 1); % length of preceding gap
            else
                ICnks{nc} = [ICnks{nc}; i]; % accumulate bins
            end
        end % of chunk,gap detecting
        %% find chunks that are long enough to save
        % keyboard
        for i=1:length(LGaps)
            ics = ICnks{i};
            lc = length(ics);
            if lc < 2 % isolated single bin
                % add this ones gap to that before next chunk
                if i<length(LGaps)
                    LGaps(i+1) = LGaps(i+1) + LGaps(i);
                end
            else
                lg = LGaps(i);
                rg = abs( ADX.depth(iG(ics(1))) - zRef ); % distance from ADCP to chunk start
                ir = find(Rng_USE(:,1)<rg & rg<=Rng_USE(:,2));
                if ~isempty(ir)
                    gspecs = Bins_USE{ir};
                    % find minimum allowed chunk size based on gap length  
                    ig = find(gspecs(:,2)<=lg & lg<=gspecs(:,3));
                    if ~isempty(ig) & lc >= gspecs(ig,1)
                        iX(ics) = 1; % flag these bins as being in good chunk
                    end
                end
                % Note that gaps for rejected chunks of 2 or more bins are
                %  not added to next gap
            end
        end % of chunk-checking (for i=1:length(LGaps))
        % Flag orphan bins (UP or DN profile)
        ix = find(~iX); % those that did not pass chunk tests
        if ~isempty(ix)
            ADX.status(iG(ix),icol) = ADX.status(iG(ix),icol) + Sta_orphan;
        end
    end % of ud=1:2 loop; UP then DN profile
end % of icol-loop; profiles
                
