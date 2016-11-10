function [latpoints2,lonpoints2]=SWIMSProfileLocs(latpoints,lonpoints,speed,dep)

%determine where each profile would be based on the simple model for SWIMS drag.
latmid=nanmean(latpoints.');
re=6.3767e+006;

	xdist=(lonpoints-lonpoints(1))*2*pi*re*cos(latmid/180*pi)/360 / 1000;
	ydist=(latpoints-latpoints(1))*2*pi*re/360 / 1000;
	dxdist=diff(xdist);
	dydist=diff(ydist);
	ds=sqrt(dxdist.^2+dydist.^2);

[dt,dx]=ComputeSWIMSspacings(speed,dep);
lonpoints2=[];
latpoints2=[];
counter=0;
for c=1:length(latpoints)-1
	%# of profiles in this leg
	numper=fix(ds(c)./dx);
    %Now make points for each profile, excluding those that actually lie on waypoints.
	latpoints2(counter+1:counter+numper-1)=latpoints(c)+(1:numper-1)/numper*(latpoints(c+1)-latpoints(c));
	lonpoints2(counter+1:counter+numper-1)=lonpoints(c)+(1:numper-1)/numper*(lonpoints(c+1)-lonpoints(c));
	counter=counter+numper-1;
end
