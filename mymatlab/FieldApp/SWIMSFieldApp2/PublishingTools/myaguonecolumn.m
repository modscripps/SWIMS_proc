function myaguonecolumn(vextent,horizoffset)
% myaguonecolumn(vextent,horizoffset)
% 
% Adjusts the paper size of the current figure to represnt the size of a 
% one column figure in AGU format.  (3 3/8 inches)
%   vextent - specifies how much of the page space the figure should take
%             vertically.  We assume a text height of 8.5 inches.  Therefore
%             vextent=1/3 implies figure height of 2.833 inches.
%	horizoffset -- default 1/2 inch.  If more is needed, figure us scooted over.
% Note that this is also a good place to set some defaults for fontsize
% and fontweights.

if ~exist('horizoffset')
	horizoffset=0.5;
end

CWIDTH = 3+3/8;
CHEIGHT = 8.5*vextent;

un=get(gcf,'units');
set(gcf,'units','inches','paperpos',[0.75 0.75 CWIDTH CHEIGHT]);
set(gcf,'units',un);

set(gcf,'defaultaxesfontsize',8);
set(gcf,'defaulttextfontsize',8);
set(gcf,'defaulttextfontsize',8);
set(gcf,'defaultaxeslinewidth',0.75);
set(gcf,'defaultlinelinewidth',1);

% set the default axes postion.  We need enough room for the
% xlabel and ylabel.  

% The axis is 0.5 inches smaller than CWIDTH 
ho=horizoffset
hs=ho+.15;
set(gcf,'defaultaxesposition',...
  [ho/CWIDTH ho/CHEIGHT (CWIDTH-hs)/CWIDTH (CHEIGHT-hs)/CHEIGHT]);

