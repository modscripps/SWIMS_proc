function mapax=PSShorePlot(latmin,latmax,lonmin,lonmax)%function mapax=PSShorePlot(latmin,latmax,lonmin,lonmax)%Plot up the PS bathymetry from Miles Logsdon - in B/W with %the shoreline heavy and 50-m contours light.[lat,lon,topo]=PSBathy;rat=(latmax-latmin)/(lonmax-lonmin);%lat/lon ratio.%indlon=find(lon > lonmin & lon < lonmax);indlat=find(lat > latmin & lat < latmax);%Antoher try -- this works great.land=50;topo(find(topo>-1))=land;[c,h]=contour(lon(indlon),lat(indlat),topo(indlat,indlon),[0 0],'k-');lw(h,2)hold on[c,h]=contour(lon(indlon),lat(indlat),topo(indlat,indlon),[-200 -150 -100 -50],'k');hold offset(gca,'XDir','reverse')%ylabel('Latitude / degrees')xlabel('West Longitude / degrees') ylabel('North Latitude / degrees') ff=.001;axis([lonmin-ff lonmax+ff latmin-ff latmax+ff ])dar=get(gca,'dataaspectratio');latmid=48;dar(1)=dar(1).*cos(latmid/180*pi);set(gca,'plotboxaspectratio',dar)set(gca,'box','on')mapax=gca;