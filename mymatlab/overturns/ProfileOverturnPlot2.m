function ax=ProfileOverturnPlot(CTDout,Overturns,PlotP);
%A summary plot of a profile.


if nargin < 3
PlotP.wh=1;
PlotP.zmin=0;
PlotP.zmax=max(CTDout.Z(:,1));
PlotP.plotcoast=1;
PlotP.plotk=0;
PlotP.plotThT=0;
end

wh=PlotP.wh;

if PlotP.zmax==0
    PlotP.zmax=max(CTDout.Z(:,wh));
end
%Allow for gridded or not by noting the size of the pressure variable
[m,n]=size(CTDout.Z);
if n==1
    whp=1;
else 
    whp=wh;
end


iz=find(CTDout.Z(:,whp)>PlotP.zmin & CTDout.Z(:,whp)<PlotP.zmax);
whs=iz(1);
whf=iz(end);
zmin=PlotP.zmin;
zmax=PlotP.zmax;

if ~isempty(Overturns)
%Find the overturns in this profile
o=ExtractStructureData(Overturns);
ind=find(o.wh_out==wh);
else
    ind=[];
end

%figure(2)
clf
ax=MySubplot(.1,.3,0,.05,.1,0.08,4,2);
axes(ax(1))

plot(CTDout.D(whs:whf,wh),CTDout.Z(whs:whf,whp),CTDout.Ds(whs:whf,wh),CTDout.Z(whs:whf,whp))
axis ij
xlabel('D')
ylim([zmin zmax])
ylabel('Depth / m')

axes(ax(2))
plot(CTDout.T(whs:whf,wh),CTDout.Z(whs:whf,whp),CTDout.Ts(whs:whf,wh),CTDout.Z(whs:whf,whp))
axis ij
ytloff
xlabel('T')
%title(['\epsilon= ' num2str(Overturn.eps) ', Lt= ' num2str(Overturn.Lt) ', GKe= ' num2str(Overturn.GKe)])
ylim([zmin zmax])
axes(ax(3))
plot(CTDout.S(whs:whf,wh),CTDout.Z(whs:whf,whp))
axis ij
ytloff
xlabel('S')
ylim([zmin zmax])
title(num2str(wh))

title(['# ' num2str(CTDout.id(wh)) ', ' num2str(CTDout.year(wh)) ', (' num2str(CTDout.lat(wh)) ', ' num2str(CTDout.lon(wh)) '), H=' num2str(CTDout.H(wh))])

axes(ax(4))
set(gca,'visible','off')
axes(ax(5))


%            plot(CTDout.Th(whs:whf,wh),CTDout.Z(whs:whf,whp),CTDout.ThT(whs:whf,wh),CTDout.Z(whs:whf,whp),...
%                CTDout.ThC(whs:whf,wh),CTDout.Z(whs:whf,whp))
tmp=abs(CTDout.ThT(whs:whf,wh));
tmp(find(tmp>100))=100;
tmp(find(tmp<.1))=.1;
%plot(abs(CTDout.Th(whs:whf,wh)),CTDout.Z(whs:whf,whp),-tmp,CTDout.Z(whs:whf,whp))
if PlotP.plotThT==1
    semilogx(abs(CTDout.Th(whs:whf,wh)),CTDout.Z(whs:whf,whp),tmp,CTDout.Z(whs:whf,whp),'m-')
else
    semilogx(abs(CTDout.Th(whs:whf,wh)),CTDout.Z(whs:whf,whp))
end
axis ij
hold on
for d=1:length(ind)
    c=ind(d);
    %h=plot([1e-10 1e-10],[Overturns(c).zs Overturns(c).zf],'b-');
    err=Overturns(c).err;
    h=plot([Overturns(c).Lt Overturns(c).Lt],[Overturns(c).zs Overturns(c).zf],'r-');
    lw(h,2)
    h2=plot([err.Lz err.Lz],[Overturns(c).zs Overturns(c).zf],'g-');
    lw(h2,2)
    h2=plot([err.Lr err.Lr],[Overturns(c).zs Overturns(c).zf],'b-');
    lw(h2,2)

end
hold off

%xlim([-100 100])
xlim([1e-1 1e2])
xlabel('Thorpe')
%ytloff
ylim([zmin zmax])

axes(ax(6))
semilogx(CTDout.N2(whs:whf,wh),CTDout.Z(whs:whf,whp))
axis ij
ytloff
hold on
for d=1:length(ind)
    c=ind(d);
    h=plot([Overturns(c).N2 Overturns(c).N2],[Overturns(c).zs Overturns(c).zf],'r-');
    lw(h,2)
end
hold off
xlabel('N^2')
%xlim([0 1.2*max(CTDout.N2(whs:whf,wh))])
xlim([1e-7 1e-3])
ylim([zmin zmax])
set(gca,'xtick',[1e-6 1e-4])
%EPSILON
%n2v=[1e-7 1e-6 1e-5 1e-4 1e-3];
axes(ax(7))
epsmax=1e-4;
tmp=CTDout.eps(:,wh);
tmp(find(tmp>epsmax))=epsmax;
%semilogx(tmp,CTDout.Z(:,whp),'k-')
h=semilogx(tmp,CTDout.Z(:,whp),'k-');
lw(h,6)

hold on
for d=1:length(ind)
    c=ind(d);
    err=Overturns(c).err;
    h=plot([Overturns(c).eps Overturns(c).eps],[Overturns(c).zs Overturns(c).zf],'r-');
    lw(h,2)
    h=plot([err.epsz err.epsz],[Overturns(c).zs Overturns(c).zf],'g-');
    lw(h,2)
    h=plot([err.epsr err.epsr],[Overturns(c).zs Overturns(c).zf],'b-');
    lw(h,2)
end
hold off
axis ij
set(gca,'xtick',[1e-10 1e-8 1e-6 1e-4])
xlim([1e-10 epsmax])
ylim([zmin zmax])
grid
ytloff
xlabel('\epsilon')
set(gca,'xscale','log')


axes(ax(8))
if PlotP.plotk==1
kmax=1e-1;
tmp=0.2*CTDout.eps(:,wh)./CTDout.N2(:,wh);
tmp(find(tmp>epsmax))=kmax;
h=semilogx(tmp,CTDout.Z(:,whp),'k-');
lw(h,4)
hold on
for d=1:length(ind)
    c=ind(d);
    err=Overturns(c).err;
    h=plot([0.2*Overturns(c).eps./Overturns(c).N2 0.2*Overturns(c).eps./Overturns(c).N2],[Overturns(c).zs Overturns(c).zf],'r-');
    lw(h,2)
    h=plot([err.Kz err.Kz],[Overturns(c).zs Overturns(c).zf],'g-');
    lw(h,2)
    h=plot([err.Kr err.Kr],[Overturns(c).zs Overturns(c).zf],'b-');
    lw(h,2)

end
hold off
axis ij
set(gca,'xtick',[1e-6 1e-4 1e-2])
xlim([1e-7 kmax])
ylim([zmin zmax])
grid
ytloff
xlabel('K')
set(gca,'xscale','log')
else %plot GKe and run len
hold on
for d=1:length(ind)
    c=ind(d);
    h=plot([Overturns(c).GKe Overturns(c).GKe],[Overturns(c).zs Overturns(c).zf],'g-');
    lw(h,2)
    h2=plot([Overturns(c).GKrunlen Overturns(c).GKrunlen],[Overturns(c).zs Overturns(c).zf],'b-');
    lw(h2,2)
end
plot([CTDout.PP.CUTOFF CTDout.PP.CUTOFF],[zmin zmax],'b--')
plot([CTDout.PP.THRESH_GK CTDout.PP.THRESH_GK],[zmin zmax],'g--')
axis ij
set(gca,'xtick',[1e-1 1e0 1e1 1e2])
xlim([1e-1 1e2])
hold off
ylim([zmin zmax])
grid
ytloff
xlabel('GKe')
set(gca,'xscale','log')
if ~isempty(ind)
lg=legend([h(1) h2(1)],'e','run');
set(lg,'fontsize',8)
end    
    
end


axes('position',[.75 .2 .23 .23])
plot(CTDout.S(whs:whf,wh),CTDout.T(whs:whf,wh))
axis square
xlabel('S')
ylabel('T')
if ~isempty(ind)
%Plot another diagnostic of the overturn scale versus n2.
%n2=g/rho*drho/dz.  So dz= g/rho/n2 * delrho.
n2v=[1e-7 1e-6 1e-5 1e-4 1e-3];
err2=GKerr(sqrt(n2v),0,err.drho,err.dz,CTDout.PP.n);
%K=0.2eps/N^2 =0.2(0.64Lt^2N^3)/N^2 so
%% that Lt=sqrt(10 K / N). K=1e-4 in deep ocean (N=5e-4 rad/s) gives Lt=1.4
%o=ExtractOverturnData(Overturns);

axes('position',[.75 .6 .20 .24])
% h=loglog(o.N2(ind),o.Lt(ind),'r.',n2v,err2.Lz,'g--',n2v,err2.Lr,'b--')
 h=loglog(o.N2(ind),o.Lt(ind),'r.',n2v,err2.Lz,'g--',n2v,err2.Lr,'b--',...
     n2v,err2.ltsig3,'m--',n2v,err2.ltsig4,'m-')
 

h=legend(h,'Data',['\Delta z=' num2str(err.dz)],['\Delta \rho=' num2str(err.drho)],...
    'K=10^{-3}','K=10^{-4}',3)
set(h,'fontsize',8)
 xlim([1e-7 1e-3])
set(gca,'xtick',[1e-7 1e-5 1e-3])
xlabel('N^2')
ylabel('L_T')
end

if PlotP.plotcoast==1
    axes('position',[.75 .85 .20 .14])
    load coast
    h=plot(long,lat,npi2pi(CTDout.lon(wh)),CTDout.lat(wh),'r.');
    set(h(end),'markersize',25)
end    
xtloff
ytloff
