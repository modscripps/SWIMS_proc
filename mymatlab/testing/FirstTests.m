%yb=70.3; ye=70.9; % First test, off shelf, soak and zero H2S on deck, then cycle
yb=72.0; ye=72.5; % N-S line on mar 14th
CV = get_SWIMS_RawData(yb, ye, 'D:\swims\BS03\indexes\CTD_BS03_matfiles.mat',...
    'D:\swims\BS03\data_mat\CTD', 1);
plot(CV.yday_adj,CV.Pr,'k-',CV.yday_adj,CV.Pr,'r.'), grid on
axis ij
title('Pr vs Yday, initial SWIMS tests in BS03')
figure
plot(CV.yday_adj,AtoD_SWIMS(CV.addata(6,:),2)), grid on
title('H2S Volts')
figure
plot(CV.yday_adj,AtoD_SWIMS(CV.addata(8,:),2)), grid on
title('pH Volts')


return
% during Grid_NewCTD.m to check start,end indices for down/up cycles
K>> plot(CSci.Pr(IND)), hold on
K>> plot(inds,CSci.Pr(IND(inds)),'g.',indf,CSci.Pr(IND(indf)),'r.')
K>> axis ij, zoom on, grid on
K>> figure
K>> plot(inds,'g.'), hold on, plot(indf,'r.'), zoom on, grid on
