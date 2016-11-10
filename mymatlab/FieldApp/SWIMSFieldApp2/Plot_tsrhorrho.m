function [axts,axrrho]=Plot_tsrhorrho(SWIMSgrid,wh,pmin,pmax,pos_tsd)
ip=find(SWIMSgrid.z_ctd>pmin & SWIMSgrid.z_ctd < pmax);

S1=SWIMSgrid.s1(ip,wh)*1000;
D1=SWIMSgrid.sgth1(ip,wh);
T1=SWIMSgrid.t1(ip,wh);

P=SWIMSgrid.z_ctd(ip);


tmin=min(T1);tmax=max([tmin+.01 max(T1)]);
smin=min(S1);smax=max([smin+.01 max(S1)]);
dmin=min(D1);dmax=max([dmin+.01 max(D1)]);

Tlimits=[tmin tmax pmin pmax];
Slimits=[smin smax pmin pmax];
Dlimits=[dmin dmax pmin pmax];
limits=[Tlimits; Slimits;Dlimits];
%
X=[T1 S1 D1];
Y=P*ones(1,3);
xlabeltext=['         T / {}^o C        ';
    '         S / psu           ';
    '\sigma_{\theta} / kg m^{-3}']; 
titletext='';
ylabeltext='p  / dbar  ';
linetype=[' r-';' g-';' b-'];
dy=.07;

axts=mymultixaxis3(X,Y,[0 0 0],xlabeltext,ylabeltext,titletext,[],[],[],[],[],...
    pos_tsd,dy,limits,1,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Now plot rrho too
[Rrho,angle]=MyRrhofun2(S1,T1,P/100,8);
pgrid=P;
%Rrho
isf=find(Rrho>1);
idr=find(Rrho>0 & Rrho < 1);
ist=find(Rrho < 0);

xmin=-2;
xmax=2;

%clip at xmin and xmax
rrp=Rrho;
rrp(find(Rrho>xmax))=xmax;
rrp(find(Rrho<xmin))=xmin;
%axes('Position',pos_rrho)
axrrho=axes;
plot(rrp,pgrid,rrp(isf),pgrid(isf),'b.',rrp(idr),pgrid(idr),'r.')
axis ij
ylim([pmin pmax])
xlim([xmin xmax])
set(gca,'XTick',(xmin:xmax))
set(gca,'YTicklabel','')
grid
xlabel('R_{\rho}')
