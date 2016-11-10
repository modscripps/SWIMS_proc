function [ah,ahcb]=ContourSWIMSdata(SWIMSgrid,varstr,beg_time,end_time,zmin,zmax,dmin,dmax,sm)
%function [ah,ahcb]=ContourSWIMSdata(SWIMSgrid,varstr,beg_time,end_time,zmin,zmax,dmin,dmax)


x=SWIMSgrid.yday;
y=SWIMSgrid.z_ctd;
%varstr='sgth1';

tind=find(SWIMSgrid.yday>beg_time & SWIMSgrid.yday<end_time);
zind=find(SWIMSgrid.z_ctd>zmin & SWIMSgrid.z_ctd<zmax);
minmaxd=[dmin dmax];

eval(['data=SWIMSgrid.' varstr '(zind,tind);'])

[ah,ahcb]=Bitchinplot2(x(tind),y(zind),data,minmaxd,'jet',varstr,2,1,2,sm,0,1,32);

axis([beg_time end_time zmin zmax])
ylabel('Depth / m')
xlabel('Yearday')