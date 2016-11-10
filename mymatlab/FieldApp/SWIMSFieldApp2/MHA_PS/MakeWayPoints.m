function MakeWayPoints(speed,dep)
%function MakeWayPoints(speed,dep)
%MakeWaypoints.m
%Interactively plot waypoints on the current figure.
%It uses the figure axis limits to determine the mid latitude.
%Compute the total distance and time elapsed.
%Also overplot the approximate spacing of profiles, given speed in knots and 
%depth of profiling in m.  These default to 4 knots and 200 m depth if not specified.
%4/02 MHA
%
if nargin==0
    speed=4;
    dep=200;
end

%IN M/S
speedms=speed/2;
%in km/hr
speedkmh=speed*1.6;

%Get the mean lat
latmid=mean(ylim);

%Earth radius
re=6.3767e6;

lonpoints=[];
latpoints=[];
np=1;
[tmpx,tmpy]=ginput(1);
while ~isempty(tmpx)
	lonpoints(np)=tmpx;
	latpoints(np)=tmpy;
	hold on
	h=plot(lonpoints,latpoints,'b');
	lw(h,2)
	h=plot(lonpoints,latpoints,'r.');
	set(h,'MarkerSize',8)
	hold off
	%sum up the total distance in km
	xdist=(lonpoints-lonpoints(1))*2*pi*re*cos(latmid/180*pi)/360 / 1000;
	ydist=(latpoints-latpoints(1))*2*pi*re/360 / 1000;
	dxdist=diff(xdist);
	dydist=diff(ydist);
	ds=sqrt(dxdist.^2+dydist.^2);
	totdist=sum(ds);
    tottime=totdist/speedkmh;
	title(['Speed: ' num2str(speed) ' kts, cycle depth: ' num2str(dep) 'm, Total distance: ' num2str(totdist) ' km, time =' num2str(tottime) ' hr'])
	[tmpx,tmpy]=ginput(1);
	np=np+1;
end
[latpoints' lonpoints']

latmins=(latpoints-floor(latpoints))*60;
lonmins=(lonpoints-floor(lonpoints))*60;
[latmins' lonmins']

[latpoints2,lonpoints2]=SWIMSProfileLocs(latpoints,lonpoints,speed,dep);
hold on
	h=plot(lonpoints2,latpoints2,'g.');
	set(h,'MarkerSize',8)
hold off
%title('')
