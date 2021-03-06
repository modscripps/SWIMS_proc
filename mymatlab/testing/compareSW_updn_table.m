% compareSW_updn_table
% e.g. vars = {'ts12','Ts12', 1 if vs. depth, 'xlabel', 'ylabel'}

vars = {'ts12','Ts12',0,'S * 1000','T / ^{\circ}C'; ...
        'tc12','Tc12',0,'C / S m^{-1}','T / ^{\circ}C'; ...
        't1dox','T1Dox',0,'Dox / ml L^{-1}','T1 / ^{\circ}C'; ...
        't2dox','T2Dox',0,'Dox / ml L^{-1}','T2 / ^{\circ}C'; ...
        
        'dt12','Dt12',1,'T_1 - T_2 / ^{\circ}C', []; ...
        'dc12','Dc12',1,'C_1 - C_2 / S m^{-1}',[]; ...
        'ds12','Ds12',1,'(S_1 - S_2) * 1000',[]; ...
        'dth12','Dth12',1,'\theta_1 - \theta_2 / ^{\circ}C',[]; ...
        'dsg12','Dsg12',1,'\sigma_{\theta}_1 - \sigma_{\theta}_2 / kg m^{-3}',[]; ...
        
        't12','T12',1,'T_1, T_2 / ^{\circ}C',[]; ...
        'c12','C12',1,'C_1, C_2 / S m^{-1}',[]; ...
        's12','S12',1,'(S_1, S_2) * 1000',[]; ...
        'th12','Th12',1,'\theta_1, \theta_2 / ^{\circ}C',[]; ...
        'sg12','Sg12',1,'\sigma_{\theta}_1, \sigma_{\theta}_2 / kg m^{-3}',[]; ...
        
        't1','T1',1,'T_1 / ^{\circ}C',[]; ...
        't2','T2',1,'T_2 / ^{\circ}C',[]; ...
        'c1','C1',1,'C_1 / S m^{-1}',[]; ...
        'c2','C2',1,'C_2 / S m^{-1}',[]; ...
        's1','S1',1,'S_1',[]; ...
        's2','S2',1,'S_2',[]; ...
        'th1','Th1',1,'\theta_1 / ^{\circ}C',[]; ...
        'th2','Th2',1,'\theta_2 / ^{\circ}C',[]; ...
        'sgth1','Sg1',1,'\sigma_{\theta}_1 / kg m^{-3}',[]; ...
        'sgth2','Sg2',1,'\sigma_{\theta}_2 / kg m^{-3}',[]; ...
        'obs','Obs',1,'Obs / FTU',[]; ...
        'obs2','Obs2',1,'Obs2 / FTU',[]; ...
        'flu','Flu',1,'Flu / {\mu}g L^{-1}',[]; ...
        'dox','Dox',1,'Dox / ml L^{-1}',[]; ...
        'h2s','H2S',1,'H2S / {\mu}g L^{-1}',[]; ...
        'pH','pH',1,'pH',[]; ...
        'sulf','Sulf',1,'Sulf / {\mu}g L^{-1}',[]};

cruz = {'ps01';'ps02'; 'ps03'; 'bs03'; 'hc03'; 'ml04'; 'aeg04'; 'stf07';
    'philex08';'mc09';'mort';'ArcticMix'};