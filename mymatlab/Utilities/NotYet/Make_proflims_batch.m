% Make_proflims_batch.m

clear
set_swims_paths
SWIMSfoldersPC

CRUs = {'ps02','ps01'};

for CRUn = 1:2
    CRU = CRUs{CRUn};
    indF = ['c:/swims/' CRU '/indexes/CTD_' CRU '_matfiles.mat'];
    matF = ['c:/swims/' CRU '/data_mat/CTD'];
    
    load(indF)
    
    clear PROFS;
    iPro = 0; ybE = NaN;
    
    for iMf = 1:length(Index.yday_beg)
        yb = Index.yday_beg(iMf) - 2/1440;
        ye = Index.yday_end(iMf) + 1/1440;
        if ~isnan(ybE) & (yb-ybE)*1440<10
            yb = ybE - 1/1440;
        end
        CTD = get_SWIMS_RawData(yb, ye, indF, matF); 
        PRO = find_Swims_profiles(CTD.Pr, CTD.yday_adj);
        if ~isempty(PRO)
            iPro = iPro+1;
            PROFS(iPro) = PRO;
            % start next set just before last profile of this one
            if ~isempty(PRO.yday_beg)
                ybE = PRO.yday_beg(end);
            end
        end
    end
    
    save(['ProfLims_' CRU], 'PROFS')
end

CRU='ps02';
load(['ProfLims_' CRU])
indF = ['c:/swims/' CRU '/indexes/CTD_' CRU '_matfiles.mat'];
matF = ['c:/swims/' CRU '/data_mat/CTD'];

for iP=1:length(PROFS)
    %figure(2-mod(iP,2)), clf
    clear CTD; iP
    CTD = get_SWIMS_RawData(PROFS(iP).yday_beg(1)-2/1440, PROFS(iP).yday_end(end)+2/1440, ...
        indF, matF);
    % Take out code thru 'return' to plot profiles, this is to find weird points
    ix = find(mod(diff(CTD.modCT),256)~=1 | ...
        diff(CTD.yday_adj)<=0 | diff(CTD.yday_adj)>(2/86400) );
    if ~isempty(ix)
        figure(1)
        subplot(3,1,1)
        plot(1:length(CTD.modCT),CTD.modCT,'k-',1:length(CTD.modCT),CTD.modCT,'g.',...
            ix,CTD.modCT(ix),'ro')
        title(['yday_beg = ' num2str(CTD.yday_adj(1))]), grid on
        CTD.yday_adj = CTD.yday_adj - floor(CTD.yday_adj(1));
        subplot(3,1,2)
        plot(1:length(CTD.modCT),CTD.yday_adj,'k-',1:length(CTD.modCT),CTD.yday_adj,'b.',...
            ix,CTD.yday_adj(ix),'ro'), grid on
        subplot(3,1,3)
        plot(1:length(CTD.modCT),CTD.Pr,'k-',1:length(CTD.modCT),CTD.Pr,'c.',...
            ix,CTD.Pr(ix),'ro'), axis ij, grid on
        zoom on
        pause
    end
end
return

    plot(CTD.yday_adj(1:4:end),CTD.Pr(1:4:end),'k-'), axis ij, hold on, zoom on
    
    plot(PROFS(iP).yday_beg,PROFS(iP).pr_beg,'g.',...
        PROFS(iP).yday_end,PROFS(iP).pr_end,'r.')
    for ia=1:length(PROFS(iP).yday_beg)
        c='y-'; l=2.5;
        if abs(PROFS(iP).type(ia))>0 & abs(PROFS(iP).pr_end(ia)-PROFS(iP).pr_beg(ia))>.2
            l=1;
        end
        switch (PROFS(iP).type(ia))
        case 1
            c='g-';
        case -1
            c='m-';
        case 0
            c='b-';
        end
        plot([PROFS(iP).yday_beg(ia) PROFS(iP).yday_end(ia)], ...
            [PROFS(iP).pr_beg(ia) PROFS(iP).pr_end(ia)], c, 'linewidth',l)
    end
    zoom on
    pause
end
