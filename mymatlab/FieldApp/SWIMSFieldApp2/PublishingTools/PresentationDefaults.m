disp('changing defaults...')
disp(' ')
%% General Figure defaults:
set(0,'DefaultFigureColor','w','DefaultFigureInvertHardCopy','off');
set(0,'DefaultSurfaceEdgeColor','k','DefaultSurfaceFaceColor','interp');
set(0,'DefaultFigurePaperOrientation','portrait');
set(0,'DefaultFigurePaperPosition',[0.25 2.6111 8 5.7778]);
set(0,'DefaultFigurePosition',[25 125 768 555]);
%% Axes, Labels, and Numbering:
%set(0,'DefaultAxesView',[-20 30],'DefaultAxesColor','w');
set(0,'DefaultAxesXColor','k','DefaultAxesYColor','k','DefaultAxesZColor','k')
if (get(0,'ScreenDepth') == 1)
  set(0,'DefaultAxesColorOrder',[0 0 0]);
else
  %set(0,'DefaultAxesColorOrder',[0 0 0; .9 0 0; 0 .5 0; 0 0 1; 0 .5 .7; .7 0 .7; .7 .5 0; .6 .5 .6]);
  set(0,'DefaultAxesColorOrder',[0 0 0; 1 0 0; 0 0 1; 0 1 0; 0 .5 .7; .7 0 .7; .7 .5 0; .6 .5 .6]);
end
set(0,'DefaultAxesFontName','Times New Roman','DefaultAxesFontSize',14)

set(0,'DefaultAxesFontWeight','Bold')

set(0,'DefaultAxesLineWidth',2)

set(0,'DefaultAxesTickLength',[.02 .04])

set(0,'DefaultAxesLineStyleOrder','-|--|:|-.')
%                             y   x  xw  yh
set(0,'DefaultAxesPosition',[0.1 0.1 0.8 0.8]);

%% Other text:
set(0,'DefaultTextColor',[0 0 0]);


