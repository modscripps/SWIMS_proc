% Run_SWIMS_sync

cruise = 'ArcticMix'; %'wa_nliw_apr2013'; %'MORT'; %'mc09'; %'ps02';

while 1
    xx = inputdlg('yd_b,yd_e', 'CTD/ADCP sync', 1,{''});
    if isempty(xx) || isempty(str2num(xx{1}))
        break
    end
    yd = str2num(xx{1});
    yd_b=yd(1); yd_e=yd(2);

    AD_CTD_Synch
end

%% You can make this a loop, where you enter new yd_b,yd_e 
%  to look at another interval