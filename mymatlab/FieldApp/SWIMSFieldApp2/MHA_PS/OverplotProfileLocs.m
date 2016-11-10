function h=OverplotProfileLocs(latpoints,lonpoints,speed,dep)
%function h=OverplotProfileLocs(latpoints,lonpoints,speed,dep)
%On the current figure, overplot the locations of SWIMS up/down profiles
%given a set of waypoints, a speed in knots and a cycle depth in m.
%Return a handle to the plot.
%4/02 MHA
%

%First compute the locations of the profiles given the speed and depth
[latpoints2,lonpoints2]=SWIMSProfileLocs(latpoints,lonpoints,speed,dep)
%Then plot them.
hold on
	h=plot(-lonpoints2,latpoints2,'g.');
	set(h,'MarkerSize',8)
hold off
%title('')
