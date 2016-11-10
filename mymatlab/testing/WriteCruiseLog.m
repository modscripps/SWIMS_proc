function WriteEventLog

%%%%%%%%%%%%%%%%% INFO TO EDIT FOR EACH CRUISE %%%%%%%%%%%%%%%%%%%%%%%%%%
year = 2003;
get_gps = 'y';

cruise.name = 'bs03';
cruise.year = year;


                                            
%%%%%%%%%%%%%%%%%%%%%%%%%%%% INFO THAT SHOULDN'T NEED EDITING %%%%%%%%%%%%%%
cruise.log.script = [cruise.name '\matlab\logs\WriteCruiseLog'];
cruise.log.MAX_TIME_DIFF = 10/(24*3600);  % Maximum yday difference allowed to 
                                            %accept gps position of drop  
YDAY_INCR = 100/(24*3600); % Increment subtracted from 1st drop time & added 
                           % to last drop time before getting gps data
cruise.log.datestr = datestr(clock);


% load <cruise>_folders
eval([cruise.name '_folders'])

cruise.leg(1).start=yearday(10,3,cruise.year,16,6);

cruise.stn(1).name = 'slope1';
cruise.stn(1).yday = [yearday(10, 3, year, 21) yearday(18, 3, year, 4, 30)];
cruise.stn(2).name = 'crossgyre1';
cruise.stn(2).yday = [yearday(18, 3, year, 9, 2) yearday(19, 3, year, 7, 21)];
cruise.stn(3).name = 'crossgyre2';
cruise.stn(3).yday = [yearday(19, 3, year, 11, 31) yearday(20, 3, year, 4, 40)];
cruise.stn(4).name = 'midgyre';
stn(4).yday = [yearday(20, 3, year, 6, 5) yearday(23, 3, year, 14)];
%cruise.stn(5).name = 'slope2';
%cruise.stn().yday = [yearday(, year,) yearday(, year,)];
n_stn = length(cruise.stn);

%%%%%%%%%%%%%%%%%%%%%%%% AMP log info %%%%%%%%%%%%%%%%%%%%%%%%
cruise.amp.burst(1).drops = [18500:18502];
cruise.amp.burst(1).waypts ='';
cruise.amp.burst(2).drops = [18503:18512];
cruise.amp.burst(2).waypts ='CrossSlope1';
cruise.amp.burst(3).drops = [18513:18518];
cruise.amp.burst(3).waypts ='CrossSlope1';
cruise.amp.burst(4).drops = [18519:18527];
cruise.amp.burst(4).waypts ='CrossSlope1';
cruise.amp.burst(5).drops = [18528:18534];
cruise.amp.burst(5).waypts ='CrossSlope1';
cruise.amp.burst(6).drops = [18535:18542];
cruise.amp.burst(6).waypts ='CrossSlope1';
cruise.amp.burst(7).drops = [18543:18550];
cruise.amp.burst(7).waypts ='CrossSlope1';
cruise.amp.burst(8).drops = [18551:18562];
cruise.amp.burst(8).waypts ='CrossSlope1';
cruise.amp.burst(9).drops = [18563:18574];
cruise.amp.burst(9).waypts ='CrossSlope1';
cruise.amp.test(1).drops = 18575;
cruise.amp.burst(10).drops = [18576:18579];
cruise.amp.burst(10).waypts ='AlongSlope1';
cruise.amp.burst(11).drops = [18580:18586];
cruise.amp.burst(12).waypts ='AlongSlope4';
cruise.amp.burst(12).drops = [18587:18590];
cruise.amp.burst(12).waypts ='AlongSlope4';
cruise.amp.burst(13).drops = [18591:18602];
cruise.amp.burst(13).waypts ='AlongSlope4';
cruise.amp.burst(14).drops = [18603:18607];
cruise.amp.burst(15).waypts ='AlongSlope5';
cruise.amp.burst(15).drops = [18608:18617];
cruise.amp.burst(15).waypts ='AlongSlope5';
cruise.amp.burst(16).drops = [18618:18623];
cruise.amp.burst(16).waypts ='AlongSlope5';
cruise.amp.burst(17).drops = [18624:18636];
cruise.amp.burst(17).waypts ='AlongSlope5';
cruise.amp.burst(18).drops = [18637:18641];
cruise.amp.burst(18).waypts ='AlongSlope5';
cruise.amp.burst(19).drops = [18642:18647];
cruise.amp.burst(19).waypts ='AlongSlope5';
cruise.amp.burst(20).drops = [18648:18653]; % start CrossGyre1 stn
cruise.amp.burst(21).drops = [18654:18659];
cruise.amp.burst(22).drops = [18660:18664];
cruise.amp.burst(23).drops = [18665:18667];
cruise.amp.burst(24).drops = [18668:18669];
cruise.amp.burst(24).comment = 'Depth limited to 500 m';
cruise.amp.burst(25).drops = [18670:18675]; % start CrossGyre2 stn
cruise.amp.burst(26).drops = [18676:18679];
cruise.amp.burst(27).drops = [18680:18684]; % start midgyre stn
cruise.amp.burst(28).drops = [18685:18691];
cruise.amp.burst(29).drops = [18692:18695];
%cruise.amp.burst().drops = [];
%cruise.amp.burst().waypts ='';

% check amp inputs for consistency %%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:length(cruise.amp.burst)           
    if  cruise.amp.burst(j).drops(1) > cruise.amp.burst(j).drops(end)
        display(['amp burst ' int2str(j) ' drops out of sequence'])
    end
end

eval(['load ' char(39) amp_procdata '\amplog' char(39)])

% Add a vector 'yday' to cruise.amp.burst(j) with the times of all drops in
% the burst.
for j = 1:length(cruise.amp.burst) 
    % First, get yday for this drop and add it to structure
    yday = []; 
    for k = 1:length(cruise.amp.burst(j).drops)
        i_drop = find(amplog(:,1) == cruise.amp.burst(j).drops(k));
        yday = [yday amplog(i_drop,3)];
    end
    cruise.amp.burst(j).yday = yday;
    
    % If get_gps is set to 'y', get position for each drop and add it to
    % structure
    lat = []; lon = [];
    if strcmp(get_gps, 'y') | strcmp(get_gps, 'yes')
        GPSdata = get_gps_data(yday(1)-YDAY_INCR, yday(end)+YDAY_INCR, ...
                        gps_indexfile, gps_datapath);
        % Calculate the minimum yday difference between the drop time and
        % sattime in the gps data.  If the difference is less than a preset
        % yday increment, take lat and lon for that index value.
        for k = 1:length(cruise.amp.burst(j).drops)
            dt = min(abs(GPSdata.sattime - cruise.amp.burst(j).yday(k)));
            if dt <= cruise.log.MAX_TIME_DIFF
                i_dt = find(abs(GPSdata.sattime - cruise.amp.burst(j).yday(k)) == dt);
                lat = [lat GPSdata.lat(i_dt(1))];
                lon = [lon GPSdata.lon(i_dt(1))];
            else
                lat = [lat NaN];
                lon = [lon NaN];
            end
        end
    end
    cruise.amp.burst(j).lat = lat;
    cruise.amp.burst(j).lon = lon;
    
    % add fields 'stn_name' and 'stn_number' to cruise.amp.burst(j)
    m = 1; mm = 1;
    while mm == 1 & m < n_stn
        if cruise.amp.burst(j).yday(1) <= amplog(i_drop,3) & ...
                cruise.amp.burst(j).yday(end) >= amplog(i_drop,3)  
            cruise.amp.burst(j).stn_name  = cruise.stn(m).name;
            cruise.amp.burst(j).stn_number = m;
            mm = 0;
        end
        m = m + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%% Box Core info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cruise.boxcore(1).yday = yearday(11,3,year,6,28);
cruise.boxcore(1).lat = 41 + 23.448/60; cruise.boxcore(1).lon = 30 + 0.624/60;

%%%%%%%%%%%%%%%%%%%%%%%% CTD info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cruise.ctd(1).yday = yearday(10,3,year,21);
cruise.ctd(1).lat = 41 + 32.338/60; cruise.ctd(1).lon = 30 + 0.459/60;
cruise.ctd(1).comment = 'To get sound speed for SeaBeam mapping';
cruise.ctd(2).yday = yearday(11,3,year,8,30);
cruise.ctd(2).lat = 41 + 32.338/60; cruise.ctd(2).lon = 30 + 00.459/60;
cruise.ctd(2).comment = 'rosette';
%cruise.ctd(1).yday
cruise.ctd(3).lat = 41 + 21.324/60; cruise.ctd(3).lon = 30 + 0.730/60;
cruise.ctd(3).comment = 'rosette';

%%%%%%%%%%%%%%%%%%%%%%%%% MMP info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stn = 'slope1';
cruise.mmp.burst(1).drops = [13900:13909];
cruise.mmp.burst(2).drops = [13910:13945];
cruise.mmp.burst(2).waypts = 'AlongShelf1';
cruise.mmp.burst(3).drops = [13946:13960];
cruise.mmp.burst(3).waypts = 'AlongShelf1';
cruise.mmp.burst(4).drops = [13961:13971]; % start midgyre stn
cruise.mmp.burst(5).drops = [13972:13983];
cruise.mmp.burst(6).drops = [13990:13996];
cruise.mmp.burst(7).drops = [13977:14012];
cruise.mmp.burst(8).drops = [14013:14026];
cruise.mmp.burst(9).drops = [14027:14034];
cruise.mmp.burst(10).drops = [14035:14044];
cruise.mmp.burst(11).drops = [14045:14059];
%cruise.mmp.burst().drops = [14060:14067];
%cruise.mmp.burst().drops = [14 14];
%cruise.mmp.burst().drops = [14 14];
%cruise.mmp.burst().waypts = '';

cruise.mmp.getbathy = [1:3];
cruise.mmp.getgps = [1:length(cruise.mmp.burst)];

%%%%%%%%%%%%%%%%%%%%%%%% Multibeam info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cruise.multibeam(1).yday = [yearday(10,3,year,23) yearday(12,3,year,6)];
cruise.multibeam(1).comment = 'Slope Site 1';

%%%%%%%%%%%%%%%%%%%%%%%% SWIMS info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cruise.swims.run(1).yday=[yearday(12,3,year,11,35) yearday(12,3,year,14,10)];
cruise.swims.run(2).yday=[yearday(14,3,year,4,49) yearday(14,3,year,6,40)];
cruise.swims.run(2).waypts = 'AlongSlope1';
cruise.swims.run(3).yday=[yearday(14,3,year,6,48) yearday(14,3,year,7,25)];
cruise.swims.run(3).waypts = 'CrossSlope1';
cruise.swims.run(4).yday=[yearday(14,3,year,7,43) yearday(14,3,year,8,39)];
cruise.swims.run(4).waypts = 'CrossSlope1';
cruise.swims.run(5).yday=[yearday(14,3,year,8,49) yearday(14,3,year,10,15)];
cruise.swims.run(5).waypts = 'CrossSlope1';
cruise.swims.run(6).yday=[yearday(14,3,year,10,15) yearday(14,3,year,11,15)];
cruise.swims.run(6).waypts = 'CrossSlope1';
cruise.swims.run(7).yday=[yearday(14,3,year,11,15) yearday(14,3,year,12,26)];
cruise.swims.run(7).waypts = 'CrossSlope1';
cruise.swims.run(8).yday=[yearday(14,3,year,12,30) yearday(14,3,year,13,45)];
cruise.swims.run(8).waypts = 'CrossSlope1';
cruise.swims.run(9).yday=[yearday(14,3,year,13,52) yearday(14,3,year,14,59)];
cruise.swims.run(9).waypts = 'CrossSlope1';
cruise.swims.run(10).yday=[yearday(14,3,year,15,02) yearday(14,3,year,16,27)];
cruise.swims.run(10).waypts = 'CrossSlope1';
cruise.swims.run(11).yday=[yearday(14,3,year,16,43) yearday(14,3,year,17,41)];
cruise.swims.run(11).waypts = 'CrossSlope1';
cruise.swims.run(12).yday=[yearday(14,3,year,17,52) yearday(14,3,year,18,54)];
cruise.swims.run(12).waypts = 'CrossSlope1';
cruise.swims.run(13).yday=[yearday(14,3,year,19,02) yearday(14,3,year,20,03)];
cruise.swims.run(13).waypts = 'CrossSlope1';
cruise.swims.run(14).yday=[yearday(14,3,year,20,08) yearday(14,3,year,21,05)];
cruise.swims.run(14).waypts = 'CrossSlope1';
cruise.swims.run(15).yday=[yearday(14,3,year,21,09) yearday(14,3,year,22,04)];
cruise.swims.run(15).waypts = 'CrossSlope1';
cruise.swims.run(16).yday=[yearday(14,3,year,22,12) yearday(14,3,year,23,12)];
cruise.swims.run(16).waypts = 'CrossSlope1';
cruise.swims.run(17).yday=[yearday(14,3,year,23,14) yearday(15,3,year,00,28)];
cruise.swims.run(17).waypts = 'CrossSlope1';
cruise.swims.run(18).yday=[yearday(15,3,year,0,39) yearday(15,3,year,01,36)];
cruise.swims.run(18).waypts = 'CrossSlope1';
cruise.swims.run(19).yday=[yearday(15,3,year,01,43) yearday(15,3,year,02,48)];
cruise.swims.run(19).waypts = 'CrossSlope1';
cruise.swims.run(20).yday=[yearday(15,3,year,02,50) yearday(15,3,year,03,47)];
cruise.swims.run(20).waypts = 'CrossSlope1';
cruise.swims.run(21).yday=[yearday(15,3,year,23,10) yearday(16,3,year,0,40)];
cruise.swims.run(21).waypts = 'AlongSlope4';
cruise.swims.run(22).yday=[yearday(16,3,year,0,50) yearday(16,3,year,2,51)];
cruise.swims.run(22).waypts = 'AlongSlope4';
cruise.swims.run(23).yday=[yearday(16,3,year,3,43) yearday(16,3,year,5,23)];
cruise.swims.run(23).waypts = 'AlongSlope5';
cruise.swims.run(24).yday=[yearday(16,3,year,5,34) yearday(16,3,year,7,1)];
cruise.swims.run(24).waypts = 'AlongSlope5'; 
cruise.swims.run(24).comment = 'Begin holding swims 100 m above bottom during run';
cruise.swims.run(25).yday=[yearday(16,3,year,7,11) yearday(16,3,year,9,9)];
cruise.swims.run(25).waypts = 'AlongSlope5';
cruise.swims.run(25).comment = 'Holding swims 100 m above bottom during run';
cruise.swims.run(26).yday=[yearday(16,3,year,9,21) yearday(16,3,year,10,42)];
cruise.swims.run(26).waypts = 'AlongSlope5';
cruise.swims.run(26).comment = {'Holding swims 100 m above bottom during run.' ...
                                  'Slow loop at 1042 to avoid traffic.'};
cruise.swims.run(27).yday=[yearday(16,3,year,10,56) yearday(16,3,year,13,16)];
cruise.swims.run(27).waypts = 'AlongSlope5';
cruise.swims.run(25).comment = 'Holding swims 100 m above bottom during run';
cruise.swims.run(28).yday=[yearday(16,3,year,13,31) yearday(16,3,year,14,38)];
cruise.swims.run(28).waypts = 'AlongSlope5';
cruise.swims.run(28).comment = 'Holding swims 100 m above bottom during run';
cruise.swims.run(29).yday=[yearday(16,3,year,14,43) yearday(16,3,year,15,50)];
cruise.swims.run(29).waypts = 'AlongSlope5';
cruise.swims.run(29).comment = 'Holding swims 100 m above bottom during run';
cruise.swims.run(29).waypts = 'AlongSlope5';
cruise.swims.run(30).yday=[yearday(16,3,year,16,2) yearday(16,3,year,16,57)];
cruise.swims.run(30).waypts = 'AlongSlope5';
cruise.swims.run(31).yday=[yearday(16,3,year,17,3) yearday(16,3,year,18,0)];
cruise.swims.run(31).waypts = 'AlongSlope5';
cruise.swims.run(32).yday=[yearday(17,3,year,8,56) yearday(17,3,year,10,1)]; % end time uncertain
cruise.swims.run(32).waypts = 'Slope1Survey1';
cruise.swims.run(32).comment = 'leg 1';
cruise.swims.run(33).yday=[yearday(17,3,year,10,6) yearday(17,3,year,11,43)];
cruise.swims.run(33).waypts = 'Slope1Survey1';
cruise.swims.run(33).comment = 'leg 2';
cruise.swims.run(34).yday=[yearday(17,3,year,11,48) yearday(18,3,year,12,31)];
cruise.swims.run(34).waypts = 'Slope1Survey1';
cruise.swims.run(34).comment = 'leg 3';
cruise.swims.run(35).yday=[yearday(17,3,year,12,36) yearday(17,3,year,13,52)];
cruise.swims.run(35).waypts = 'Slope1Survey1';
cruise.swims.run(35).comment = 'leg 4';
cruise.swims.run(36).yday=[yearday(17,3,year,13,58) yearday(17,3,year,14,41)];
cruise.swims.run(36).waypts = 'Slope1Survey1';
cruise.swims.run(36).comment = 'leg 5';
cruise.swims.run(37).yday=[yearday(17,3,year,14,43) yearday(17,3,year,15,52)];
cruise.swims.run(37).waypts = 'Slope1Survey';
cruise.swims.run(37).comment = 'leg 6';
cruise.swims.run(38).yday=[yearday(17,3,year,16,2) yearday(17,3,year,16,56)];
cruise.swims.run(38).waypts = 'Slope1Survey2';
cruise.swims.run(38).comment = 'leg 1';
cruise.swims.run(39).yday=[yearday(17,3,year,17,1) yearday(17,3,year,17,24)];
cruise.swims.run(39).waypts = 'Slope1Survey2';
cruise.swims.run(39).comment = 'leg 2';
cruise.swims.run(40).yday=[yearday(17,3,year,17,25) yearday(17,3,year,18,21)];
cruise.swims.run(40).waypts = 'Slope1Survey2';
cruise.swims.run(40).comment = 'leg 3';
cruise.swims.run(41).yday=[yearday(17,3,year,18,25) yearday(17,3,year,18,43)];
cruise.swims.run(41).waypts = 'Slope1Survey2';
cruise.swims.run(41).comment = 'leg 4';
cruise.swims.run(42).yday=[yearday(17,3,year,18,44) yearday(17,3,year,19,42)];
cruise.swims.run(42).waypts = 'Slope1Survey2';
cruise.swims.run(42).comment = 'leg 5';
cruise.swims.run(43).yday=[yearday(18,3,year,10,40) yearday(18,3,year,12,50)];
cruise.swims.run(43).waypts = 'Star1';
cruise.swims.run(43).comment = 'leg 1';
cruise.swims.run(44).yday=[yearday(18,3,year,12,54) yearday(18,3,year,15,14)];
cruise.swims.run(44).waypts = 'Star1';
cruise.swims.run(44).comment = 'leg 2';
cruise.swims.run(45).yday=[yearday(18,3,year,15,22) yearday(18,3,year,17,38)];
cruise.swims.run(45).waypts = 'Star1';
cruise.swims.run(45).comment = 'leg 3';
cruise.swims.run(46).yday=[yearday(18,3,year,17,41) yearday(18,3,year,19,44)];
cruise.swims.run(46).waypts = 'Star1';
cruise.swims.run(46).comment = 'leg 4';
cruise.swims.run(47).yday=[yearday(18,3,year,19,48) yearday(18,3,year,22,2)];
cruise.swims.run(47).waypts = 'Star1';
cruise.swims.run(47).comment = 'leg 5';
cruise.swims.run(48).yday=[yearday(19,3,year,11,31) yearday(19,3,year,12,56)];
cruise.swims.run(48).waypts = 'Star2';
cruise.swims.run(48).comment = 'leg 1';
cruise.swims.run(49).yday=[yearday(19,3,year,13,2) yearday(19,3,year,14,24)];
cruise.swims.run(49).waypts = 'Star2';
cruise.swims.run(49).comment = 'leg 2';
cruise.swims.run(50).yday=[yearday(19,3,year,14,29) yearday(19,3,year,15,58)];
cruise.swims.run(50).waypts = 'Star2';
cruise.swims.run(50).comment = 'leg 3';
cruise.swims.run(51).yday=[yearday(19,3,year,16,4) yearday(19,3,year,17,26)];
cruise.swims.run(51).waypts = 'Star2';
cruise.swims.run(51).comment = 'leg 4';
cruise.swims.run(52).yday=[yearday(19,3,year,17,31) yearday(19,3,year,18,54)];
cruise.swims.run(52).waypts = 'Star2';
cruise.swims.run(52).comment = 'leg 5';
cruise.swims.run(53).yday=[yearday(20,3,year,6,19) yearday(20,3,year,7,37)];
cruise.swims.run(53).waypts = 'Star3';
cruise.swims.run(53).comment = 'leg 1';
cruise.swims.run(54).yday=[yearday(20,3,year,7,44) yearday(20,3,year,8,5)];
cruise.swims.run(54).waypts = 'Star3';
cruise.swims.run(54).comment = 'leg 2, interrupted by bad signal in tensionmeter';
cruise.swims.run(55).yday=[yearday(20,3,year,9,35) yearday(20,3,year,11,12)];
cruise.swims.run(55).waypts = 'Star3';
cruise.swims.run(55).comment = 'leg 2';
cruise.swims.run(56).yday=[yearday(20,3,year,11,18) yearday(20,3,year,12,40)];
cruise.swims.run(56).waypts = 'Star3';
cruise.swims.run(56).comment = 'leg 3';
cruise.swims.run(57).yday=[yearday(20,3,year,12,46) yearday(20,3,year,13,17)];
cruise.swims.run(57).waypts = 'Star3';
cruise.swims.run(57).comment = 'leg 4';
cruise.swims.run(58).yday=[yearday(20,3,year,14,14) yearday(20,3,year,14,45)];
cruise.swims.run(58).waypts = 'Star3';
cruise.swims.run(58).comment = 'leg 5';
cruise.swims.run(59).yday=[yearday(23,3,year,0,14) yearday(23,3,year,1,12)];
cruise.swims.run(59).waypts = 'MidGyreGrid';
cruise.swims.run(59).comment = 'leg 1';
cruise.swims.run(6).yday=[yearday(23,3,year,) yearday(23,3,year,)];
cruise.swims.run(6).waypts = 'MidGyreGrid';
cruise.swims.run(6).comment = 'leg ';
cruise.swims.run(6).yday=[yearday(23,3,year,) yearday(23,3,year,)];
cruise.swims.run(6).waypts = 'MidGyreGrid';
cruise.swims.run(6).comment = 'leg ';
cruise.swims.run(6).yday=[yearday(23,3,year,) yearday(23,3,year,)];
cruise.swims.run(6).waypts = 'MidGyreGrid';
cruise.swims.run(6).comment = 'leg ';
cruise.swims.run(6).yday=[yearday(23,3,year,) yearday(23,3,year,)];
cruise.swims.run(6).waypts = 'MidGyreGrid';
cruise.swims.run(6).comment = 'leg ';
cruise.swims.run(6).yday=[yearday(23,3,year,) yearday(23,3,year,)];
cruise.swims.run(6).waypts = 'MidGyreGrid';
cruise.swims.run(6).comment = 'leg ';
cruise.swims.run(6).yday=[yearday(23,3,year,) yearday(23,3,year,)];
cruise.swims.run(6).waypts = 'MidGyreGrid';
cruise.swims.run(6).comment = 'leg ';
cruise.swims.getbathy = [1:42];
cruise.swims.getgps = [1:length(cruise.swims.run)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XCP info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cruise.xcp(1).drop = 202;
cruise.xcp(1).yday = yearday(20,3,year,7,24);
cruise.xcp(1).lat = 42 + 30.369/60;
cruise.xcp(1).lon = 30 + 32.720/60;
cruise.xcp(2).drop = 203;
cruise.xcp(2).yday = yearday(20,3,year,13,26);
cruise.xcp(2).lat = 42 + 29.878/60;
cruise.xcp(2).lon = 30 + 31.288/60;
cruise.xcp(3).drop = 204;
cruise.xcp(3).yday = yearday(20,3,year,19,14);
cruise.xcp(3).lat = 42 + 31.146/60;
cruise.xcp(3).lon = 30 + 31.146/60;
cruise.xcp(4).drop = 205;
cruise.xcp(4).yday = yearday(21,3,year,1,20,30);
cruise.xcp(4).lat = 42 + 30.201/60;
cruise.xcp(4).lon = 30 + 30.265/60;
cruise.xcp(5).drop = 206;
cruise.xcp(5).yday = yearday(21,3,year,7,30,43);
cruise.xcp(5).lat = 42 + 31.346/60;
cruise.xcp(5).lon = 30 + 29.442/60;
cruise.xcp(6).drop = 208;
cruise.xcp(6).yday = yearday(21,3,year,13,9,54);
cruise.xcp(6).lat = 42 + 27/60;
cruise.xcp(6).lon = 30 + 30/60;
cruise.xcp(7).drop = 209;
cruise.xcp(7).yday = yearday(21,3,year,19,30,28);
cruise.xcp(7).lat = 42 + 29.368/60;
cruise.xcp(7).lon = 30 + 30.271;;
cruise.xcp(8).drop = 210;
cruise.xcp(8).yday = yearday(22,3,year,1,26);
cruise.xcp(8).lat = 42 + 31.470/60;
cruise.xcp(8).lon = 30 + 29.764/60;
cruise.xcp(9).drop = 211;
cruise.xcp(9).yday = yearday(22,3,year,7,26,0);
cruise.xcp(9).lat = 42 + 29.002/60;;
cruise.xcp(9).lon = 30 + 29.998/60;
cruise.xcp(10).drop = 212;
cruise.xcp(10).yday = yearday(22,3,year,14,3);
cruise.xcp(10).lat = 42 + 27.843/60;
cruise.xcp(10).lon = 30 + 30.224/60;
cruise.xcp(11).drop = 213;
cruise.xcp(11).yday = yearday(22,3,year,19,43,39),;
cruise.xcp(11).lat = 42 + 32.087/60;
cruise.xcp(11).lon = 30 + 29.534/60;
cruise.xcp(12).drop = 214;
cruise.xcp(12).yday = yearday(23,3,year,1,20,5);
cruise.xcp(12).lat = 42 + 32.533/60;
cruise.xcp(12).lon = 30 + 25.897/60;
%cruise.xcp().drop = 21;
%cruise.xcp().yday = yearday(21,3,year,);
%cruise.xcp().lat = 42 + ;
%cruise.xcp().lon = 30 + ;
%cruise.xcp().drop = 21;
%cruise.xcp().yday = yearday(21,3,year,);
%cruise.xcp().lat = 42 + ;
%cruise.xcp().lon = 30 + ;
%cruise.xcp().drop = 21;
%cruise.xcp().yday = yearday(21,3,year,);
%cruise.xcp().lat = ;
%cruise.xcp().lon = ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK INPUTS FOR CONSISTENCY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

chk_mmp = 'y';
chk_swims = 'y';



 
%%%%%%%%%%%%%%%%%%%%%%%%% check mmp %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(chk_mmp, 'y') 
    for j = 1:length(cruise.mmp.burst)           
        if  cruise.mmp.burst(j).drops(1) > cruise.mmp.burst(j).drops(end)
                display(['mmp burst ' int2str(j) ' drops out of sequence'])
         end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%% check swims %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(chk_swims, 'y') 
    for j = 1:length(cruise.swims.run)           
        if  cruise.swims.run(j).yday(1) > cruise.swims.run(j).yday(end)
            display(['swims run ' int2str(j) ' ydays are out of sequence'])
        end
    end
end



eval(['save ' char(39) 'C:\' cruise.name '\matlab\logs\' cruise.name '_cruiselog' char(39) ' cruise'])