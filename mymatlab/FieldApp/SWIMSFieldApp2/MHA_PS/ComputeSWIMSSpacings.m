function [dt,dx]=ComputeSWIMSSpacings(speed,dep)
%[dt,dx]=ComputeSWIMSSpacings(speed,dep)
%Specify speed in knots and maximum desired depth for SWIMS to cycle
%dt will be the time in seconds between profiles, up and down counted as separate profiles.
%dx will be the distance in km.
%4/02 MHA


%This is a matrix of ship speed in m/s and the resultant sine of the angle - 
%so that this is depth/line out.  Note this is computed for a constant 
%line out, 600 m, so it is only an approximation.
ShipSinAngle=   [[0.2500    0.9955
    0.5000    0.9456
    0.7500    0.8486
    1.0000    0.7533
    1.2500    0.6747
    1.5000    0.6115
    1.7500    0.5602
    2.0000    0.5181
    2.2500    0.4827
    2.5000    0.4527
    2.7500    0.4270
    3.0000    0.4044
    3.2500    0.3847
    3.5000    0.3672
    3.7500    0.3515
    4.0000    0.3373
    4.2500    0.3242
    4.5000    0.3127
    4.7500    0.3021
    5.0000    0.2919]];

%IN M/S
speedms=speed/2;
%in km/hr
speedkmh=speed*1.6;


%Now get the (approx) line out needed
val=interp1(ShipSinAngle(:,1),ShipSinAngle(:,2),speedms);
lineout=dep/val;

%Mean fall rate, averaged between down (3/2 m/s) and up (1 m/s).
fallrate=1.25;
%time between profiles, up and down counted each as profiles
dt=lineout/fallrate;
%distance between profiles, in km
dx=speedms*dt/1000;
