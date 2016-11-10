function OverturnCloseupPlot(CTDout,Overturn)
%A simple summary plot of an overturn.
dzb=10;

whs=Overturn.s-dzb;
whf=Overturn.f+dzb;
if whs<1
    whs=1;
end
if whf>length(CTDout.D)
    whf=length(CTDout.D);
end

wh=Overturn.wh_out;
[m,n]=size(CTDout.Z);
if n==1
    whp=1;
else 
    whp=wh;
end


%figure(2)
clf
ax=MySubplot(.1,.3,0,.1,.1,0,4,1);
axes(ax(1))

plot(CTDout.D(whs:whf,wh),CTDout.Z(whs:whf,whp),CTDout.Ds(whs:whf,wh),CTDout.Z(whs:whf,whp))
axis ij
xlabel('D')

axes(ax(2))
plot(CTDout.T(whs:whf,wh),CTDout.Z(whs:whf,whp),CTDout.Ts(whs:whf,wh),CTDout.Z(whs:whf,whp))
axis ij
ytloff
xlabel('T')
title(['\epsilon= ' num2str(Overturn.eps) ', Lt= ' num2str(Overturn.Lt) ', GKe= ' num2str(Overturn.GKe)])
axes(ax(3))
plot(CTDout.S(whs:whf,wh),CTDout.Z(whs:whf,whp))
axis ij
ytloff
xlabel('S')

axes(ax(4))
plot(CTDout.Th(whs:whf,wh),CTDout.Z(whs:whf,whp),CTDout.ThT(whs:whf,wh),CTDout.Z(whs:whf,whp),...
    CTDout.ThC(whs:whf,wh),CTDout.Z(whs:whf,whp))
axis ij
xlabel('Thorpe')
axes('position',[.72 .4 .24 .24])
plot(CTDout.S(whs:whf,wh),CTDout.T(whs:whf,wh),CTDout.S(whs+dzb:whf-dzb,wh),CTDout.T(whs+dzb:whf-dzb,wh))
axis square
