% Examine_RDI_VmDas.m  -  convert raw RDI file to matlab, compute useful quantities,
%    display various parameters, contour profile data

clear Vel RD dd
set_swims_paths
disp('Select raw RDI file to examine:')
[FilNm,PthNm] = uigetfile([pwd '\ML*','Raw RDI file');

if ~FilNm
    disp('None selected, exiting')
    return
end
ok=1;
if length(FilNm)<5
    ok=0;
elseif strcmp(FilNm(end-3:end), '.mat')
    ok=0;
end
if ~ok
    disp('Invalid file selected, start over.')
    return
end

Rf = fullfile(PthNm,FilNm);
OutDir = tempdir;
if isempty(OutDir)
    OutDir = pwd;
end
Mf = fullfile(OutDir,['Mat' FilNm(3:end)]);
disp([' Matlab file = ' Mf])

raw2mat_V3(Rf, Mf, [], {'.mat'}); % To look at all RDI data (mainly header)
Vel = Get_ADCP_SwimsDN(Rf, -1, 2); % To look at computed quantities
RD=load(Mf);

% Gather fixed header values and display
clear subs % organize into subtypes
subs(1).typ = '-SYSTEM-';
subs(1).flds = {'sysconfig','firmwareversion','firmwarerevison',...
        'numberbeams','sensors','sensorson'};
subs(2).typ = '-MISC-';
subs(2).flds = {'coords','headoffset','headbias','lowcorr','goodthresh',...
        'errthresh','fishthresh','watermode','waterband'};
subs(3).typ = '-SAMPLING-';
subs(3).flds = {'npings','nbins','binlen','blanklen','dis1',...
        'pulselen','pulselag','codereps'};
subs(4).typ = '-BOTTOM Track-';
if isempty(RD.btmaxdep)
    subs(4).flds = {'OFF'};
else
    subs(4).flds = {'btnpings','btmaxdep','btmode','btdelay',...
            'btminamp','btcorrthresh','btmaxerr'};
end
cc = {' FIELD ','   Value'};
for i=1:length(subs)
    cc(end+1,1:2) = { subs(i).typ, '========' };
    for ia=1:length(subs(i).flds)
        if isfield(RD, subs(i).flds{ia})
            x = RD.(subs(i).flds{ia});
            if length(x)>1, x=x(1); end
        else
            x = ' missing';
        end
        cc(end+1,1:2) = { subs(i).flds{ia}, x };
    end
end
disp(cc)

fE = figure;
plot(Vel.yday,Vel.ens_no,'k.');
title('Yearday vs. Ensemble no.');
fA = figure;
plot(Vel.yday,Vel.pitch,'g-', Vel.yday,Vel.roll,'r-'); title('pitch(g),roll(r)');
figure(fE)
disp(['Check Settings listing above, and examine Figure ' int2str(fE) ' and ' int2str(fA) ';'])
disp('Press any key when ready for further plots: ')
pause
fB = figure;
figure(fA), clf
plot(Vel.yday, Vel.SWIMS_headT,'k-'); title('Heading (T)');
xx = int2str(fA);
if ~isempty(RD.btmaxdep)
    figure(fB)
    plot(Vel.yday,Vel.VbtN,'g', Vel.yday,Vel.VbtE,'r'); title('BTvel E(r) N(g)');
    xx = [xx ' and ' int2str(fB)];
end

x=input(['Examine Figure ' xx '; Then press <enter> to start contour plots:'],'s');
clear PW
PW(1) = figure(fA); clf;PW(2) = figure(fB); clf;PW(3) = figure; clf
disp(['For each beam, examine Figures ' int2str(PW) ' .'])
disp('When finished, press any key to check next beam.')
for ib=1:4
    disp(['Beam ' int2str(ib) ' contours ... '])
    eval(['dd=Vel.v' num2str(ib) '_bm;'])
    figure(PW(1)),clf,pcolor(Vel.yday,Vel.z_adcp,dd), hold on
    shading flat,colormap(jet),colorbar,axis ij
    caxis([-2.5 2.5]), colorbar
    axis([Vel.yday(1) Vel.yday(end) 0 110])
    xlabel('yearday'),ylabel('depth / m')
    title(['Vel' num2str(ib) ' (beam,m/s)'])
    plot(Vel.yday, Vel.bottomBT,'k-','linewidth',2.5)
    ;
    eval(['dd=Vel.ec' num2str(ib) '_bm;'])
    figure(PW(2)),clf,pcolor(Vel.yday,Vel.z_adcp,dd), hold on
    shading flat,colormap(jet),colorbar,axis ij
    axis([Vel.yday(1) Vel.yday(end) 0 110])
    caxis([40 230]), colorbar
    xlabel('yearday'),ylabel('depth / m')
    title(['Echo' num2str(ib)])
    plot(Vel.yday, Vel.bottomBT,'k-','linewidth',2.5)
    
    eval(['dd=Vel.ec' num2str(ib) '_bm;'])
    figure(PW(3)),clf,pcolor(Vel.yday,Vel.z_adcp,dd), hold on
    shading flat,colormap(jet),colorbar,axis ij
    axis([Vel.yday(1) Vel.yday(end) 0 110])
    caxis([40 160]), colorbar
    xlabel('yearday'),ylabel('depth / m')
    title(['Cor' num2str(ib)])
    plot(Vel.yday, Vel.bottomBT,'k-','linewidth',2.5)
    
    pause
end