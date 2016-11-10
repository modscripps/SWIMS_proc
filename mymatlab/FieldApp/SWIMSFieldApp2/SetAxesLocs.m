%figure position
%fig_pos= [30 111 1306 841];
%fig_pos= [30 100 841 800];
SW.fig_pos=[35 72 768 875];
SW.fig_pos2=[ 829    64   568   430];
SW.fig_pos=[35 90 770 590];
SW.fig_pos2=[ 500    90   550   430];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%axes positions
%lat/lon plot
SW.x1=.09;
SW.y1=.7;
SW.dx1=.2;
SW.dy1=.25;

%margin
SW.ix=.03;

%colorbar
SW.x2=SW.x1+SW.dx1+SW.ix+.07;
SW.y2=SW.y1;
SW.dx2=.01;
SW.dy2=SW.dy1;

%tide plots
SW.x3=SW.x1;
SW.y3=.58;
SW.dx3=SW.dx1+SW.dx2+SW.ix+.07;
SW.dy3=.07;

SW.x4=SW.x1;
SW.y4=.5;
SW.dx4=SW.dx3;
SW.dy4=SW.dy3;

%tsd and rrho windows
SW.pos_tsd=[SW.x1 .23 .15 .18];
SW.pos_rrho=[.25 .23 .13 .18];

%contour window 1
SW.nx=3;
SW.x5=SW.x2+SW.dx2+SW.nx*SW.ix;
SW.y5=.7;
SW.dx5=.38;
SW.dy5=.25;

SW.x6=SW.x5+SW.dx5+SW.nx*SW.ix;
SW.y6=SW.y5;
SW.dx6=.01;
SW.dy6=SW.dy5;

%contour window 2
SW.x7=SW.x2+SW.dx2+SW.nx*SW.ix;
SW.y7=.4;
SW.dx7=SW.dx5;
SW.dy7=.25;

SW.x8=SW.x7+SW.dx7+SW.nx*SW.ix;
SW.y8=SW.y7;
SW.dx8=.01;
SW.dy8=SW.dy7;

%contour window 3
SW.x9=SW.x2+SW.dx2+SW.nx*SW.ix;
SW.y9=.1;
SW.dx9=SW.dx5;
SW.dy9=.25;

SW.x10=SW.x9+SW.dx9+SW.nx*SW.ix;
SW.y10=SW.y9;
SW.dx10=.01;
SW.dy10=SW.dy9;


