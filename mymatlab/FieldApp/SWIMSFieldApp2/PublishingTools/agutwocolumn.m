function aguonecolumn(vextent)
% aguonecolumn(vextent)
% 
% Adjusts the paper size of the current figure to represnt the size of a 
% one column figure in AGU format.  (3 3/8 inches)
%   vextent - specifies how much of the page space the figure should take
%             vertically.  We assume a text height of 8.5 inches.  Therefore
%             vextent=1/3 implies figure height of 2.833 inches.
%
% Note that this is also a good place to set some defaults for fontsize
% and fontweights.

CWIDTH = 7;
CHEIGHT = 8.5*vextent;

un=get(gcf,'units');
set(gcf,'units','inches','paperpos',[0.75 0.75 CWIDTH CHEIGHT]);
set(gcf,'units',un);

set(gcf,'defaultaxesfontsize',9);
set(gcf,'defaulttextfontsize',9);
set(gcf,'defaulttextfontsize',9);
set(gcf,'defaultaxeslinewidth',0.75);
set(gcf,'defaultlinelinewidth',1);

% set the default axes postion.  We need enough room for the
% xlabel and ylabel.  

% The axis is 0.5 inches smaller than CWIDTH 
set(gcf,'units','norm');
set(gcf,'defaultaxesposition',...
  [0.5/CWIDTH 0.5/CHEIGHT (CWIDTH-0.65)/CWIDTH (CHEIGHT-0.65)/CHEIGHT]);

