timetowait = 180; % in s

while (true)
	
	% ?
    NewCtdFLG = 1;
	
    while (NewCtdFLG==1)
        Dout = GetNewDataSIO_func('SproulTest', [] ,year);
        NewCtdFLG = Dout.NewCtdFLG;
    end

    Dout = Grid_RunNewProfs_func('SproulTest');
% %     make_SWIMS_plots
%     pause(0)
% % %     timetorun = now+timetowait/3600/24;
% %     disp([datestr(now) ': Waiting until ' datestr(timetorun)])
% %     pause(timetowait)
end
