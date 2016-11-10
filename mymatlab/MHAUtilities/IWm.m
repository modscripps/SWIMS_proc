function m=IWm(f,k,N,lat)
%function m=IWm(f,k,N,lat)
%Internal wave vertical wavenumber in cpm, given desired frequency f
%in cph (CYCLES per hour), horizontal wavenumber k in cpm, Vaisala
%frequency N in cph.
%latitude lat in degrees.  Default lat is 30 degrees.  Default N is 3 cph.

if nargin < 2
	disp 'USAGE: IWm=m(f,k,N,lat)'
	return
elseif nargin == 2
	lat=30;
	N=3;
elseif nargin == 3
	lat=30;
elseif nargin > 4
	disp 'USAGE: IWm=m(f,k,N,lat)'
	return
end

fo=2* (1/24)*sin(lat*pi/180);

m=k.*sqrt( (N.^2 - f.^2) ./ (f.^2 - fo.^2) );

