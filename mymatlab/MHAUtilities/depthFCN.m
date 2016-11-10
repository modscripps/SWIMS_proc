%function y=depth(p,lat)

function y=depth(p,lat)
%Return depth in meters given pressure in decibars and latitude in degrees
%pressure may be a vector or a matrix; lat must be one number
%Adapted from UNESCO.c

	x=sin(lat/57.29578);
	x=x.*x;      
	gr=9.780318*(1.0+(5.2788E-3 + 2.36E-5*x)*x) + 1.092E-6*p;
    d=(((-1.82E-15*p+2.279E-10).*p -2.2512E-5).*p + 9.72659).*p ; 
	y=d./gr;


%void
%DEPTH(float p,float lat,float *depth)
%{
%/*********************************************************** 
%        UNITS:
%         PRESSURE      P    DECIBARS
%         LATITUDE      LAT  DEGREES 
%         DEPTH         DEPTH  METERS
%  
%        CHECKVALUE: DEPTH=9712.653 FOR P=10000, LAT=30
%***********************************************************/ 
%      extended		LAT,P,X,GR,D;
%	  P=p;
%	  LAT=lat;
%      X= sin(LAT/57.29578);
%      X= X*X;
%      GR=9.780318*(1.0+(5.2788E-3 + 2.36E-5*X)*X) + 1.092E-6*P;
%      D=(((-1.82E-15*P+2.279E-10)*P -2.2512E-5)*P + 9.72659)*P ; 
%      *depth=D/GR;
%	  return;
%}
