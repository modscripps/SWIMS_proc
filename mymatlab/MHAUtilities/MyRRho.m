wh=20;
str=['Da Sheet:Ri Thorpe ALL:Number ' num2str(wh)];
%cd to next directory
eval(['cd ''',str,''''])

RecordsToDo=180;
SBINS=4;
if wh > 0
	StartRec=(wh-1)*180+28;
else
	StartRec=(wh+2)*180 + 90 + 4;
end

MBL_Info;

load T
S=loadpart('Da Wattah:Final Mat Files:CTD.mat',[2+2*2816+1 2+3*2816 StartRec+1 StartRec+180]);
load D
Tm=mean(T')';
Sm=mean(S')';
Dm=mean(D')';

Tt=T(:,100);
St=S(:,100);
Dt=D(:,100);

clear T S D
load Z
Zm=mean(Z')';
Zt=Z(:,100);
clear Z

[b,a]=butter(4,1/2/13);
Zmf=filtfilt(b,a,Zm);
Tmf=filtfilt(b,a,Tm);
Smf=filtfilt(b,a,Sm);
Dmf=filtfilt(b,a,Dm);

s=700;
f=2700;

alph=-1./1026.*diff(Dmf)./diff(Tmf);
beta=1./1026.*diff(Dmf)./diff(Smf);

Rp=mean(alph(s:f))./mean(beta(s:f)).*diff(Tmf)./diff(Smf);

subplot(1,3,1)
plot(Tm(30:2816),Zm(30:2816))
axis ij
xlabel('T (¡C)')

subplot(1,3,2)
plot(Sm(30:2816),Zm(30:2816))
axis ij
xlabel('S (ppt)')

subplot(1,3,3)
plot(Dm(30:2816),Zm(30:2816))
axis ij
xlabel('D (kgm^{-3}-1000)')

suptitle('MBL 12 hour mean T,S,D profiles')

sd=1500;fd=f;
subplot(1,3,1)
plot(Tt(sd:fd),Zt(sd:fd))
axis ij
xlabel('T (¡C)')

subplot(1,3,2)
plot(St(sd:fd),Zt(sd:fd))
axis ij
xlabel('S (ppt)')

subplot(1,3,3)
plot(Dt(sd:fd),Zt(sd:fd))
axis ij
xlabel('D (kgm^{-3}-1000)')

suptitle('MBL Typical  T,S,D profiles')

semilogy(Zmf(s:f),-Rp(s:f))
grid
xlabel('Depth(m)')
ylabel('Log_{10} R_{\rho}')
title('12 Hour Mean Density Ratio')
DateStamp

plot(Zmf(s:f),alph(s:f))
title('alpha')
xlabel('Depth(m)')
ylabel('alpha (¡C^{-1})')

semilogy(Zmf(s:f),beta(s:f))
title('beta')
xlabel('Depth(m)')
ylabel('beta (ppt^{-1})')
