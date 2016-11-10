function PlotFigureGUI(SW)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make the plot track window
%Make the figure
figure(2)
clf
set(gcf,'position',SW.fig_pos2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Experiment with using the mean temperature over a depth range
var=SW.SWIMSgrid.yday;
indz=find(SW.SWIMSgrid.z_ctd>SW.zmint&SW.SWIMSgrid.z_ctd<SW.zmaxt);
tind=find(SW.SWIMSgrid.yday>SW.beg_time & SW.SWIMSgrid.yday < SW.end_time);
eval(['var=nanmean(SW.SWIMSgrid.' SW.varstr1 '(indz,tind));'])
varmin=SW.dmin1;
varmax=SW.dmax1;
if isempty(var)
    var=NaN*ones(size(tind));
end

latv=SW.SWIMSgrid.lat(tind);
lonv=SW.SWIMSgrid.lon(tind);
varstr=['<' SW.varstr1 '> from ' num2str(SW.zmint) '-' num2str(SW.zmaxt) ' m'];
%Now plot the track in color
%function [ax1,ax2]=LatLonColorPlotFcn(latv,lonv,var,varmin,varmax,varstr,latmin,latmax,lonmin,lonmax)
[ax1,ax2]=LatLonColorPlotFcn(latv,lonv,var,varmin,varmax,varstr,SW.latmin,SW.latmax,SW.lonmin,SW.lonmax);
rp=get(ax2,'position');
rp(1)=rp(1)+.03;
set(ax2,'position',rp)
axes(ax1)
imax=length(tind);
global latpoints;
global lonpoints;
hold on
if ~isempty(latpoints) & ~isempty(lonpoints) & length(latpoints)==length(lonpoints)
h=plot(abs(lonpoints),latpoints,'k-',abs(lonpoints),latpoints,'ko');
latpoints
lonpoints
lw(h(1),2)
lc(h(1),.6*[1 1 1])
end
h=plot(-lonv(imax),latv(imax),'ko',lonv(imax),latv(imax),'kx')
hold off

title(['Alongtrack Mid-depth ' SW.varstr1])


%Make the figure
figure(1)
clf
set(gcf,'position',SW.fig_pos)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Use the time as the color
varmin=SW.beg_time;
varmax=SW.end_time;
tind=find(SW.SWIMSgrid.yday>SW.beg_time & SW.SWIMSgrid.yday < SW.end_time);

var=SW.SWIMSgrid.yday(tind);

latv=SW.SWIMSgrid.lat(tind);
lonv=SW.SWIMSgrid.lon(tind);
varstr='yearday';
%Now plot the track in color
%function [ax1,ax2]=LatLonColorPlotFcn(latv,lonv,var,varmin,varmax,varstr,latmin,latmax,lonmin,lonmax)
[ax1,ax2]=LatLonColorPlotFcn(latv,lonv,var,varmin,varmax,varstr,SW.latmin,SW.latmax,SW.lonmin,SW.lonmax);
axes(ax1)
title('Track')
%Next: tide plot
[ax3,ax4]=PS02TidePlot(SW.beg_tideplot,SW.end_tideplot,SW.beg_time,SW.end_time);

%contour window 1

[ax5,ax6]=ContourSWIMSdata(SW.SWIMSgrid,SW.varstr1,SW.beg_time,SW.end_time,SW.zmin,SW.zmax,SW.dmin1,SW.dmax1,SW.sm);

%Plot bottom data
axes(ax5)
hold on
h=plot(SW.ADCP.yday,SW.ADCP.bottomBT);lw(h,2)
hold off
xlabel('')
xtloff
%contour window 2

[ax7,ax8]=ContourSWIMSdata(SW.SWIMSgrid,SW.varstr2,SW.beg_time,SW.end_time,SW.zmin,SW.zmax,SW.dmin2,SW.dmax2,SW.sm);

%Plot bottom data
axes(ax7)
hold on
h=plot(SW.ADCP.yday,SW.ADCP.bottomBT);lw(h,2)
hold off
xlabel('')
xtloff

%contour window 3- can plot vel if desired
if strcmp(SW.varstr3(1:2),'u_')==1 | strcmp(SW.varstr3(1:2),'v_')==1
    
    [ax9,ax10]=ContourADCPdata(SW.ADCP,SW.varstr3,SW.beg_time,SW.end_time,SW.zmin,SW.zmax,SW.dmin3,SW.dmax3,SW.sm);
    
    %Plot bottom data
    axes(ax9)
    hold on
    h=plot(SW.ADCP.yday,SW.ADCP.bottomBT);lw(h,2)
    hold off
    
else
    [ax9,ax10]=ContourSWIMSdata(SW.SWIMSgrid,SW.varstr3,SW.beg_time,SW.end_time,SW.zmin,SW.zmax,SW.dmin3,SW.dmax3,SW.sm);
    
    %Plot bottom data
    axes(ax9)
    hold on
    h=plot(SW.ADCP.yday,SW.ADCP.bottomBT);lw(h,2)
    hold off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TS rho and rrho plots
wh=max(find(SW.SWIMSgrid.updown==1 & SW.SWIMSgrid.yday <= SW.end_time));
%wh=max(find(SW.SWIMSgrid.updown==1));
deepest=max(find(~isnan(SW.SWIMSgrid.t1(:,wh))));

%Try an earlier one - if this is too shallow
if SW.SWIMSgrid.z_ctd(deepest)<=30
    wh=wh-2;
    if wh<1
        wh=1;
    end
end

%If that one is too shallow too, then don't plot.
deepest=max(find(~isnan(SW.SWIMSgrid.t1(:,wh))));
if SW.SWIMSgrid.z_ctd(deepest)>=30
[axts,axrrho]=Plot_tsrhorrho(SW.SWIMSgrid,wh,SW.zmin,SW.zmax,SW.pos_tsd);
set(axrrho,'position',SW.pos_rrho)
end

set(ax1,'position',[SW.x1 SW.y1 SW.dx1 SW.dy1])
set(ax2,'position',[SW.x2 SW.y2 SW.dx2 SW.dy2])

set(ax3,'position',[SW.x3 SW.y3 SW.dx3 SW.dy3])
set(ax4,'position',[SW.x4 SW.y4 SW.dx4 SW.dy4])
set(ax5,'position',[SW.x5 SW.y5 SW.dx5 SW.dy5])
set(ax6,'position',[SW.x6 SW.y6 SW.dx6 SW.dy6])


set(ax7,'position',[SW.x7 SW.y7 SW.dx7 SW.dy7])
set(ax8,'position',[SW.x8 SW.y8 SW.dx8 SW.dy8])


set(ax9,'position',[SW.x9 SW.y9 SW.dx9 SW.dy9])
set(ax10,'position',[SW.x10 SW.y10 SW.dx10 SW.dy10])



orient tall
