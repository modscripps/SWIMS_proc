function h=sub_axis(hin,pos_units);
%
% SUB_AXIS puts a new axis on the plot in old plot data units.
%
% h=subaxis(hin,[x y dx dy]);  Where hin is the axis you want to put the new
% axis on top of, and [x y dx dy] is the position  of the left hand corner
% of your axis, and size respectively in data units of the old plot.  Useful
% for embeded plots (like plots of the tide over data)
%


xmin=pos_units(1); ymin=pos_units(2); 
deltax=pos_units(3); deltay=pos_units(4);

if deltay<=0
  error(sprintf('deltay %f < 0\n',deltay));
end;
if deltax<=0
  error(sprintf('deltax %f < 0\n',deltax));
end;

xlim=get(hin,'xlim');
ylim=get(hin,'ylim');
pos=get(hin,'position');

dirx=get(hin,'xdir');
if strcmp(dirx,'reverse')
  x = pos(1)+ pos(3)*(xlim(2)-xmin)./(xlim(2)-xlim(1));
  dx= pos(3)*(deltax./(xlim(2)-xlim(1)));
else
  x = pos(1)+ pos(3)*(xmin-xlim(1))./(xlim(2)-xlim(1));
  dx= pos(3)*(deltax./(xlim(2)-xlim(1)));
end;


diry=get(hin,'ydir');
if strcmp(diry,'reverse')
  y = pos(2)+pos(4)*(ylim(2)-ymin)./(ylim(2)-ylim(1));
  dy= pos(4)*(deltay./(ylim(2)-ylim(1)));
else
  y = pos(2)+pos(4)*(ymin-ylim(1))./(ylim(2)-ylim(1));
  dy= pos(4)*(deltay./(ylim(2)-ylim(1)));
end;


h=axes('position',[x y dx dy]);

