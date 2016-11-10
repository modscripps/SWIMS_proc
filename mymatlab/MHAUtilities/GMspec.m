function shearGM=GMspec(k,N)
%function shearGM=GMspec(k,N)
%return GM shear vertical wavenumber spectrum for vector k and
%mean stratification N in rads-1.

E=6.3e-5;
b=1300;
No=.00524;
j=3;
%N=13*2*pi/3600;

%convert beta star to cpm
ko=pi*j/b*N/No / 2 / pi;


dispGM=E*b^3/pi^2/j*(No/N)^2*1./(1+k.^2/ko^2);

highk=find(abs(k) >.1);
dispGM(highk)=E*b^3/pi^2/j*(No/N)^2*1./(1+k(highk).^2/ko^2) .* .1 ./ k(highk);
strGM=(2*pi.*k).^2 .* dispGM;

%Use GM76 -- ie eqn a2.
uGM=3*E*b^3*No^2/j/pi/2*1./(1+k.^2/ko^2);
uGM(highk)=3*E*b^3*No^2/j/pi/2*1./(1+k(highk).^2/ko^2).* .1 ./ k(highk);
shearGM=(2*pi.*k).^2 .* uGM;
%multiply by 2*pi for plotting versus cpm not rpm, and divide by (pi/2)
%to account for screwup (see Gregg and Kunze)
shearGM=shearGM*2*pi/pi*2;
