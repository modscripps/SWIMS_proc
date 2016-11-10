function compareSW_updn(SWydays, dopause, CalPars)
% generalized version of compareSW_tstc.m
% input example: compareSW_updn([2003, 70.0, 71.0], 'y')

compareSW_updn_table % table

% inputs using dialogs
[iv,ok] = listdlg('ListString',vars(:,1), 'Name','Var Selection', 'PromptString','Select one or more var(s):', 'ListSize', [160 360]);
[ic,ok] = listdlg('ListString',cruz, 'Name','Cruise Selection', 'SelectionMode','single', 'PromptString','Select a cruise:');

crz = cruz{ic, 1};
var = vars(iv, 1);

% default inputs
if nargin < 2 | isempty(dopause)
    dopause = 'y';
end
if nargin < 3
    CalPars = [];
end

set_swims_paths % edit swimsfolders for desired cruise, data locations
% if strcmp(crz, 'ml04')
%     swimsgridded = 'C:\swims_04\ml04\griddata';
%     swimsindex = 'C:\swims_04\ml04\indexes';
%     swimsmatdata = 'C:\swims_04\ml04\data_mat';
% end

yr = SWydays(1); ydb=SWydays(2); yde=SWydays(3);
[CTD] = get_swims_data(ydb, yde, fullfile(swimsindex,['SWIMS_' crz '_gridfiles.mat']), swimsgridded);

% new figures
clear fnums
for i = 1:length(iv)
    fnums(i) = figure; 
end

% time vs. depth for all gridded profiles in this interval
fZ = figure;
ylabel('depth / m'); xlabel('yearday')
box on, hold on
figure(fZ); plot(CTD.yd_byz,CTD.z_ctd),axis ij
    
% plot each profile, paired with following one
ii = 1;
while ii<length(CTD.updown)
    ud1=CTD.updown(ii); ud2=CTD.updown(ii+1);
    % determine sequence: up,down; down,up; or same
    if ud1<0
        ccs1 = 'rm';
        ttl = ['DN:'];
    else
        ccs1 = 'gc';
        ttl = ['UP:'];
    end
    ttl = [ttl sprintf('%3.4f',CTD.yday(ii)) '; '];
    if ud2<0
        ccs2 = 'rm';
        ttl = [ttl 'dn:'];
    else
        ccs2 = 'gc';
        ttl = [ttl 'up:'];
    end
    ttl = [ttl sprintf('%3.4f',CTD.yday(ii+1))];
    % if both are up or down, use dashed line for 2nd profiles
    if ud1==ud2
        cct = '--';
    else
        cct = '-';
    end
    % gather data to plot, cell arrays={1st, 2nd profile}; separate primary,secondary data
    % get 24-Hz data (up/down based on gridded limits)
    yb1 = CTD.yday_LDbeg(ii); ye1 = CTD.yday_LDend(ii);
    yb2 = CTD.yday_LDbeg(ii+1)-3/86400; ye2 = CTD.yday_LDend(ii+1);
    if (yb2-ye1)*86400<10
        yb2 = ye1; % remove short turnaround gap
    end
    Ch1 = get_swims_scidata(yb1, ye1, fullfile(swimsindex,['CTD_' crz '_matfiles.mat']), ...
        fullfile(swimsmatdata, 'CTD'), [], CalPars);
    Ch2 = get_swims_scidata(yb2, ye2, fullfile(swimsindex,['CTD_' crz '_matfiles.mat']), ...
        fullfile(swimsmatdata, 'CTD'), [], CalPars);
    zdat = {Ch1.Pr*100, Ch2.Pr*100};
    
    % gather data in cell arrays datax and datay
    for ivn = 1:length(iv)        
        clear datax datay        
        datay = zdat; % default datay
        ylab = 'depth / m';
        ydir = vars{iv(ivn),3};
        if ~isempty(vars{iv(ivn),5})
            ylab = vars{iv(ivn),5};
        end
        
        switch (vars{iv(ivn),2})
            case 'Ts12'
                datay = {Ch1.T1, Ch2.T1, Ch1.T2, Ch2.T2};
                datax = {Ch1.S1*1000, Ch2.S1*1000, Ch1.S2*1000, Ch2.S2*1000}; 
            case 'Tc12'
                datay = {Ch1.T1, Ch2.T1, Ch1.T2, Ch2.T2};
                datax = {Ch1.C1, Ch2.C1, Ch1.C2, Ch2.C2};
            case 'T1Dox'
                datay = {Ch1.T1, Ch2.T1};
                datax = {Ch1.Dox, Ch2.Dox}; 
            case 'T2Dox'
                datay = {Ch1.T2, Ch2.T2};
                datax = {Ch1.Dox, Ch2.Dox}; 
                
            case 'Dt12'
                datax = {Ch1.T1 - Ch1.T2, Ch2.T1 - Ch2.T2};
            case 'Dc12'
                datax = {Ch1.C1 - Ch1.C2, Ch2.C1 - Ch2.C2};
            case 'Ds12'
                datax = {Ch1.S1*1000 - Ch1.S2*1000, Ch2.S1*1000 - Ch2.S2*1000};
            case 'Dth12'
                datax = {Ch1.Th1 - Ch1.Th2, Ch2.Th1 - Ch2.Th2};
            case 'Dsg12'
                datax = {Ch1.Sg1 - Ch1.Sg2, Ch2.Sg1 - Ch2.Sg2};
                
            case 'T12'
                datax = {Ch1.T1, Ch2.T1, Ch1.T2, Ch2.T2};
                datay = [zdat, zdat];
            case 'C12'
                datax = {Ch1.C1, Ch2.C1, Ch1.C2, Ch2.C2};
                datay = [zdat, zdat];
            case 'S12'
                datax = {Ch1.S1*1000, Ch2.S1*1000, Ch1.S2*1000, Ch2.S2*1000};
                datay = [zdat, zdat];
            case 'Th12'
                datax = {Ch1.Th1, Ch2.Th1, Ch1.Th2, Ch2.Th2};
                datay = [zdat, zdat];
            case 'Sg12'
                datax = {Ch1.Sg1, Ch2.Sg1, Ch1.Sg2, Ch2.Sg2};
                datay = [zdat, zdat];
                
            otherwise
                if ~isfield(Ch1, vars{iv(ivn),2})
                    datax = {NaN*datay{1}, NaN*datay{2}}
                else 
                    datax = {Ch1.(vars{iv(ivn),2}), Ch2.(vars{iv(ivn),2})};
                end
        end
           
        % plot data
        figure(fnums(ivn))
        clf
        ylabel(ylab)
        xlabel(vars{iv(ivn),4})
        title(ttl)
        box on
        hold on
        h{ivn,1} = plot(datax{1}, datay{1}, [ccs1(1) '-']);  
        h{ivn,2} = plot(datax{2}, datay{2}, [ccs2(1) cct]);
        if length(datax) > 2 % for plotting T12, C12, etc.
            h{ivn,3} = plot(datax{3}, datay{3}, [ccs1(2) '-']);
            h{ivn,4} = plot(datax{4}, datay{4}, [ccs2(2) cct]);
        end 
        if vars{iv(ivn),3} % if not "vs" plot reverse the y-axis
            set(gca, 'ydir', 'reverse')
        else
            set(gca, 'ydir', 'normal')
        end
%         minx = min([datax{1} datax{2}]);
%         maxx = max([datax{1} datax{2}]);
%         xlim([minx-0.5 maxx+0.5]);        
%         miny = min([datay{1} datay{2}]);
%         maxy = max([datay{1} datay{2}]);
%         ylim([miny-5 maxy+5])
    end
    if strcmp(dopause,'y')
        pause  
        for ifn = 1:length(iv)
            for ihn = 1:size(h,2)
                delete(h{ifn,ihn})
            end
        end
    else
        ii = ii + 1; % skip next one (already plotted)
    end  
    ii = ii + 1; % for 'while'-loop
end
title(['Current SWIMS = ' num2str(ydb) ':' num2str(yde)])               