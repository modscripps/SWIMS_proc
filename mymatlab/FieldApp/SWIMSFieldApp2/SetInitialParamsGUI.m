function SW=SetInitialParamsGUI;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Lat ranges and times, pick for each experiment
%Hood canal
SW.latmin=47.6;
SW.latmax=48;
SW.lonmin=122.5;
SW.lonmax=123;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Admiralty Inlet - for testing
SW.latmin=47.8;
SW.latmax=48;
SW.lonmin=122.4;
SW.lonmax=122.7;

SW.beg_time=118.7;

SW.beg_timereq=119.5;
SW.end_timereq=SW.beg_time+2;

SW.end_time=SW.end_timereq;

SW.beg_tideplot=119;
SW.end_tideplot=120.5;

% %Intrusions
SW.latmin=47.77;
SW.latmax=47.88;
SW.lonmin=122.35;
SW.lonmax=122.52;
% 
SW.beg_time=126.738;
% 
SW.beg_timereq=126.738;
SW.end_timereq=SW.beg_time+2.5;
% 
SW.end_time=SW.end_timereq;
% 
SW.beg_tideplot=125.33;
SW.end_tideplot=129.33;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set paths and other constants
%base paths. 
SW.remotebasepath='D:\swims\ps02\';
SW.localbasepath='D:\SWIMS_MHA\ps02\';

%remote paths
SW.remotedata_pathSWIMS=[SW.remotebasepath 'griddata\'];
SW.remotedata_pathADCP=[SW.remotebasepath 'data_mat\ADCP\'];
SW.remoteindexdir=[SW.remotebasepath 'indexes\'];

%local paths
SW.data_pathSWIMS=[SW.localbasepath 'griddata'];
SW.data_pathADCP=[SW.localbasepath 'data_mat\ADCP\'];
SW.localindexdir=[SW.localbasepath 'indexes\'];

%Change this later to the desired folder for print eps files
SW.prpath=[SW.localbasepath 'figs\'];

SW.index_fileSWIMS='SWIMS_ps02_gridfiles.mat';
SW.index_fileADCP='ADCP_ps02_matfiles.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot parameters

%smoothing for color contours
SW.sm=0;

%PLOT limits
SW.zmin=0;
SW.zmax=200;

SW.dmin1=8.2;
SW.dmax1=9.0;
SW.dmin2=22;
SW.dmax2=23.4;
SW.dmin3=-.5;
SW.dmax3=.5;
SW.varstr1='t1';
SW.varstr2='sgth1';
SW.varstr3='v_wat';

%depth lims for the temperature mean in the fig 2 track plot
SW.zmint=40;
SW.zmaxt=80;

%Printing filename and incremental suffix
SW.printfilenum=1;
SW.printfilename='IntrusionsPlot';

%Now set up the axes positions in the figure window
SetAxesLocs;

