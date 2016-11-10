function pstides=LoadPSTides02(beg_time);
%Load in tides for Apr/May 2002 intrusions and LS experiments
%Structure:
% 1x2 struct array with fields:
%     Height
%     Current
%     Segment
%     Channel
%     jday_year
%     jday_PST
%     yearday_UTC
%     LatLon
%     StartDatePST

%if no arguments assume intrusions
if nargin < 1
    beg_time=125;
end

%DAVE: please modify this path to the correct one.
pstidepath='D:\SWIMS_MHA\PStide_data';

intrusionsyday=124;
if beg_time < intrusionsyday
    tidefile=fullfile(pstidepath,'ps02_tides.mat');
    wh=1; %hood canal
else
    tidefile=fullfile(pstidepath,'ps02b_tides.mat');
    wh=3; %pt wells
end

eval(['load(''' tidefile ''')'])
pstides=TIDE_OUT(wh);

