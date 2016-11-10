function [ax1,ax2]=PS02TidePlot(beg_tideplot,end_tideplot,beg_time,end_time)
%Plot a two-panel plot of the tide for Apr/May 02, with the overall plot limits
%beg_tideplot and end_tideplot.  The top panel is current and the bottom panel
%is displacement in km.  The time bracketed by beg_time and end_time is plotted 
%in heavy red.

%beg_tideplot=118;
%end_tideplot=119;

%Load in the tides
pstides=loadpstides02(beg_time);
%Compute displacement
dts=mean(diff(pstides.yearday_UTC))*24*3600;
pstides.disp=cumsum(pstides.Current)*dts;

tind=find(pstides.yearday_UTC>beg_time & pstides.yearday_UTC < end_time);
tindpl=find(pstides.yearday_UTC>beg_tideplot & pstides.yearday_UTC < end_tideplot);

ax=mysubplot(.1,.1,0,.1,.1,0.01,1,2);
axes(ax(1))
h=plot(pstides.yearday_UTC(tindpl),pstides.Current(tindpl),...
    pstides.yearday_UTC(tind),pstides.Current(tind),'r-');
if length(h)>=2
lw(h(2),2)
end
ylabel('U / ms^{-1}')
xtloff
grid on
title('Tides')

axes(ax(2))
h=plot(pstides.yearday_UTC(tindpl),pstides.disp(tindpl)/1000,...
    pstides.yearday_UTC(tind),pstides.disp(tind)/1000,'r-');
if length(h)>=2
lw(h(2),2)
end
ylabel('Disp / km')
xlabel('yearday 2002')
dmin=min(pstides.disp(tindpl)/1000)-.3;
dmax=max(pstides.disp(tindpl)/1000)+.3;
ylim([dmin dmax])
%zaptick('t')
grid on
ax1=ax(1);
ax2=ax(2);
