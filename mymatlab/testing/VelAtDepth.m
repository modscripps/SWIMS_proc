

%yd_b=86.578; yd_e=86.618;
yd_b=88.62; yd_e=88.8
SP=get_SWIMS_SciData(yd_b, yd_e, 'D:\swims\BS03\indexes\CTD_BS03_matfiles.mat', ...
    'D:\swims\BS03\data_mat\CTD', 1);
LD = get_WinchLine_data(yd_b, yd_e, 'D:\swims\BS03\indexes\LD_BS03_matfiles.mat', ...
    'D:\swims\BS03\data_mat\LD');

clear SW
SW.Pr = SP.Pr;
SW.yday = SP.yday_adj;
SW.Pitch = SP.Pitch;
SW.Temp = SP.T1;
SW.Sal = SP.S1;
clear SP
% load D:\swims\BS03\data_mat\ADDN\ADDN-2003-086-1333-42.mat
% load D:\swims\BS03\data_mat\ADUP\ADUP-2003-086-1333-42.mat
load D:\swims\BS03\data_mat\ADDN\ADDN-2003-088-1423-45.mat %ADDN-2003-088-1808-45.mat
load D:\swims\BS03\data_mat\ADUP\ADUP-2003-088-1423-45.mat %ADUP-2003-088-1808-45.mat

ptsC = 24; % number of 24-Hz points for center difn fall rate
Wpr = 100*(SW.Pr(ptsC+1:end)-SW.Pr(1:end-ptsC)) / (ptsC/24);
Wdy = (SW.yday(ptsC+1:end)+SW.yday(1:end-ptsC)) / 2;

Tm0 = Wdy(1); clf
%DelAtoS = .000292; % yday offset, ADDN back to CTD
DelAtoS = 0.000123; % yday offset, ADDN back to CTD

plot((ADDN.yday-Tm0-DelAtoS)*86400,-ADDN.VbtZ,'-',(ADDN.yday-Tm0-DelAtoS)*86400,-ADDN.VbtZ,'g.'), hold on
grid on, zoom on
plot((Wdy-Tm0)*86400,Wpr,'y-',(Wdy-Tm0)*86400,Wpr,'r.')

figure
plot(ADDN.ens_no,ADDN.SWIMS_pitch,'r-',ADUP.ens_no,ADUP.SWIMS_pitch,'g-')