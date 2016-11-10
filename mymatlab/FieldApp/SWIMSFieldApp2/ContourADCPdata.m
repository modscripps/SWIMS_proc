function [ah,ahcb]=ContouADCPdata(ADCP,varstr,beg_time,end_time,zmin,zmax,dmin,dmax,sm)
%function [ah,ahcb]=ContourSWIMSdata(ADCP,varstr,beg_time,end_time,zmin,zmax,dmin,dmax)


x=ADCP.yday;
y=ADCP.z_adcp;
%varstr='sgth1';


tind=find(ADCP.yday>beg_time & ADCP.yday<end_time);
zind=find(ADCP.z_adcp>zmin & ADCP.z_adcp<zmax);
minmaxd=[dmin dmax];

eval(['data=ADCP.' varstr '(zind,tind);'])

%now zero the data that is further out than the bottom
[m,n]=size(data);
for c=1:n
    indz=find(ADCP.z_adcp(zind)>ADCP.bottomBT(tind(c)));
    data(indz,c)=NaN;
end

%Fix up the label string for vel
if varstr(1)=='u'
    varstr='u_{wat}';
elseif varstr(1)=='v'
    varstr='v_{wat}';
end

[ah,ahcb]=Bitchinplot2(x(tind),y(zind),data,minmaxd,'jet',varstr,2,1,2,sm,0,1,32);

axis([beg_time end_time zmin zmax])
ylabel('Depth / m')
xlabel('Yearday')