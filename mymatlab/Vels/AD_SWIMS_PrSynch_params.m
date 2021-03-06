% AD_SWIMS_PrSynch_params.m  -- Set up synchronization parameters between
%       SWIMS CTD and ADDN, ADUP.  Use values found by close examination
%       via AD_CTD_Synch.m, etc.
% Supplies clock offsets DelDNtoBT, DelDNtoUP, DelDNtoSW, (DelUPtoSW opt.)
%  ensemble_number offset EnsDNtoUP, and magnetic declination MagDcl
%  to AD_SWIMS_PrSynch, given cruise (ID) and yearday range yd_b,yd_e .

%% (Need to have indications that ADDN,ADUP not synched, or ADUP only)
% offset between ADUP and ADDN ensembles (UP-DN): [first yearday, ens_no offset; ...]
EnsDNtoUP = [0, 0]; % change at indicated yearday (NaN if not synch'd) 
% Offset to add to ADDN to match SWIMS CTD: [ADDN yday, offset[sec]; ...]
DelDNtoSW = [0, 0; 366, 0]; % interp1 offset between yeardays
DelUPtoSW = []; % Non-empty only if ADUP not synched to ADDN (or ADDN missing)

switch lower(cruise)
    case 'ps02' % sequence UpW, DnBT, DnW
        DelDNtoBT = -0.25;
        DelDNtoUP = -0.67;
        MagDcl = 19.5 ; % for Puget Sound, Hood Canal, approx
        EnsDNtoUP = [0, 1;  115.49, 0;  117.5, 1];
        DelDNtoSW = ...
            [114.714, -1.7;  114.84, -2.4;  114.942, -3.0;
            115.49, -1.9;  116.01, -4.7;
            117.51, 0.7;  118.1, -1.9;
            118.54, 0.0;  119.16, -3.0;
            119.88, 0.0;  120.22, -1.7];
        % above are for HC, AI:  Main Basin PS later - dpw 8/22/08
    case 'home02'
        DelDNtoBT = -0.43;
        DelDNtoUP = 0.40; % SW = 08000
        MagDcl = 10.5 ; % for Oahu, approx
        DelDNtoSW =  ...
            [245.61, -1.0;  248.025, -15.6; % For entire Kauai Ridge period
            248.03, 0;  265.68, 0; % gap in SWIMS use (some tests)
            265.6838  -1.6;  267.0, -8.6;  267.8473, -14.6;
            267.8474, -5.6;  269.6715, -15.9 ]; % For entire Mamala Bay work
    case 'bs03'
        if yd_b>76.4 && yd_b<82.0 % period with no BTing
            DelDNtoBT = -0.1; % doesn't matter, not used
            DelDNtoUP = 0.51;
        else % while BTing
            DelDNtoBT = -0.43;
            DelDNtoUP = 0.40;
        end
%         if yd_b>93.5  % no synch, ADUP,ADDN ran independently
%             EnsDNtoUP = NaN; %% OR, map from ens_UPs to ens_DNs ?
%         end
        MagDcl = 3.5; % for Black Sea, approx
        DelDNtoSW = [93.4 -20.6;93.5 -20.6]; %BS03, last synched
    case 'hc03'
        DelDNtoBT = -0.44;
        DelDNtoUP = 0.47; % SW = 09000
        MagDcl = 19.5 ; % PS,HC approx
        EnsDNtoUP = [0, 0;  63.49, 1;  65.87, 0];
        DelDNtoSW = ...
            [298.025, -0.5;  298.33, -2.5;  298.56, -3.7;
            299.08, -6.3;  299.35, -7.5;  299.625, -9.0;
            301.025, -5.0;  301.045, -5.5;
            301.1, -6.0;  301.35, -7.5;  301.605, -8.7;
            302.11, -7.7;  302.38, -9.0;  302.45, -9.0;  302.68, -10.5];
    case 'ml04'
        if yd_b>59 && yd_e<91
            DelDNtoBT = -0.1; % doesn't matter, not used
            DelDNtoUP = 0.61; % With BT off in deep water
            MagDcl = 12.3; % for ML04 near 150W, 30N (geomag value)
            EnsDNtoUP = [0, 0;  63.49, 1;  65.87, 0]; % off by 1, ??
            DelDNtoSW = ...
                [62.045, -4.5;  62.562, -7.5; % 30N, 150.5W->149.5W
                62.75, -8.5;  63.39, -12.3; % 150W, 29.6N->30.4N
                63.49, -12.5;  64.0, -15.0;  65.0, -21.0;  65.79, -25.5; % patterns
                65.87, -13.3;  66.3, -15.7;  66.79, -18.4; % PC clock adj
                66.87, -6;  66.93, -6; % ADCP clocks adj, WB1 test
                66.945, -6.5;  67.23, -8.0;  68.11, -12.5; % back to WB0, then rough seas
                68.68, -4.2;  69.0, -6.0;  69.5, -8.6;  69.9, -11.0; % in until winch problems
                70.065, -11.7;  70.2374, -13.0]; % end of SWIMS2, long live SWIMS1 !
        else
            DelDNtoBT = -0.44;
            DelDNtoUP = 0.47; % SW = 09000
            MagDcl = 19.5 ; % Lk Wash
        end
    case 'stf07'
        DelDNtoBT = -0.1; % doesn't matter, not used
        DelDNtoUP = 0.61; % With BT off in deep water
        MagDcl = 11.5; % for stf07 near 156.75W, 30.7N (NobelTec value)
        EnsDNtoUP = [0, 0;  183.5, 113;  185.19, 0; 187.12, 2; ...
            192.305, 0; 193.21, 167; 196.76, 0; 208.17, 949; 209.23, 0];
        DelDNtoSW = ...
            [183.858, -3.5;  183.96, -3.6; % tests en route
            186.07, -6.5; 186.97, -10.0; %
            187.12, -1.5; 188.0, -4.8; 188.58, -7.0;
            188.653, -2.0; 189.63, -5.0; 190.90, -10.0;
            190.901, -1.5; 191.84, -5.0; 192.308, -7.0; 193.08, -10.0;
            193.21, -2.6; 194.2, -6.0; 195.0, -9.0; 195.8, -12.0; % 194.904 on, ADs stuck
            196.765, -2.5; 198.1, -7.2; 199.52, -13.0;
            199.95, -2.8; 201.19, -7.3; 202.20, -11.0; 203.0, -14.3; 203.71, -16.5;
            204.31, -2.7; 205.19, -6.0; 205.83, -8.5;  205.84, -2.3; 206.265, -4.0;
            206.69, -5.4; 207.42, -7.7; 207.75, -9.7;
            208.18, -10.8; 209.17, -14.2; 209.24, -14.5; 209.64, -15.8] ;
    case 'philex08'
        if yd_e<44 % BT off
            DelDNtoBT = 0; % doesn't matter, not used
            DelDNtoUP = 0.61; % With BT off in deep water
            MagDcl = -1.2; % for philex off S.Ch Sea, N. Mindoro
        elseif yd_b>=44 && yd_e<999 % BT on
            DelDNtoBT = -0.44;
            DelDNtoUP = 0.47; % SW = 09000
            MagDcl = -1.0; % for philex off Mindoro, mid,S.
            if yd_b>=57
                MagDcl = -0.8; % off mid-Panay,W
            end
        end
        EnsDNtoUP = [0, 0;  39.1, 168; 42.9, 1; 44.24, 0; ...
            55.4, 1; 57.26, 18];
        DelDNtoSW = ...
            [39.15, -3.5;  39.45, -4.5; % S China Sea
            42.92, -18.0; 43.86, -21.0; % approaching + off Mindoro
            44.24, -4.3; 44.965, -7.0; % Along Mindoro (Altim mis-set)
            46.37, 0; 46.80, -1.8; % N of Apo Reef
            46.87, 0.5; 47.14, -0.5; % Along Mindoro (Altim good)
            48.42, -5.8; 48.89, -7.5; % S,W of Apo Reef (T1,C1 clog last leg)
            49.13, -8.3; 49.33, -8.8; % Sill TW, end w/winch problems
            51.57, -4.4; 52.0, -5.8; 52.575, -8.0; % Sill TW, after personnel txfer
            55.4, -2.2; 55.56, -2.5; % NW of Panay 
            57.26, -9.3; 57.4,-10.0 ]; % W of Panay, aborted when snagged
    case 'mc09'
            DelDNtoBT = -0.44;
            DelDNtoUP = 0.47; % SW = 09000
            MagDcl = 14.1;
        EnsDNtoUP = [0,287; 96,0; 98.4,80; 99.49,0];
        DelDNtoSW = ...
            [94.77,-2.5; 95.81,-6; 96.17,-6.8; 97.22,-10.3; 98.35,-14.0; % grps 2,3,4
            98.4,-2.0; 99.46,-5.7; % grp 5, NOTE: ADDN raw times had 45 min added (GetNewByType.m)
            99.5,-5.9; 100.57,-9.5; 100.61,-12.0; 101.701,-15.6; % grps 6,7
            103.604,0.5; 104,-0.5; 104.705,2.5; 105.505,-0.2; % grps 8,9
            105.87,-2.1; 106.63,-4.4; 106.788,0.5; 107.837,-3;  % grps 10,11
            109.363,0.3; 110.425,-3.3; 110.449,-3.4; 111.08,-5.7; 111.695, -7.6; % grp 12,13,14
            111.74,-10.4; 112.32,-12.2; 112.33,-12.2; 113.528,-16.0; 114.52,-19.3; % grp 15,16,17
            114.673,-21.8; 114.796,-22.1; 114.915,-22.6; 115.158,-23.5; 115.325,1.2; 115.66,0.3]; % grp 18
    case 'mort'
        DelDNtoBT = -0.44;
        DelDNtoUP = 0.47; % SW = 09000
        MagDcl = 2.3; % NobelTec value
        EnsDNtoUP = [306.85,0; 312.19045,NaN; 314.68,0]; % usually synched, ADDN failing before wire shorted, re-term'd,...
        DelDNtoSW = ...
            [306.85,-1.3; 306.9,-1.3; % Short tests en route (Low, hanging block not working out)
            308.608,-6.6; 309.5125,-9.4; 309.5127,4.0; 310.1,2.0; % Grp2, then Grid1, chg PC clock pre-gp5,sub5
            311.0,-0.9; 311.88,-3.6; 312.1904,-5.0; 312.281,-5.2; % ADDN failed just before and during Eward run of E-W line, then CTD
            314.68,4.0; 315.2,2.5; 315.64,1.2; 316.1,-0.5; 316.54,-2.0; % Re-term'd, reset PC clock
            317.12,-4.0; 317.7,-6.0; 318.341,-8.0; 318.715,-9.5; 318.79999,-9.6; % continued clock drift(s)
            318.8,6.6; 318.925,6.0]; % after PC clock reset, for last group 25 and recovery (caused CTD file 'SWIMSFasttime' problem)
        DelUPtoSW = [312.19029,-48.0; 312.281,-48.5]; % ADUP okay while ADDN failing, then cable shorted
    case 'wa_nliw_apr2013'
        DelDNtoBT = -0.44;
        DelDNtoUP = 0.47; % SW = 09000
        MagDcl = 19.0; % JM chart value
        EnsDNtoUP = [113.9,0];
        DelDNtoSW = ...
            [113.9,-1.5; 114.209,-2.2; % start of X-canyon,  before power outage (recover missing hr later)
            114.26,7.9; 115.04,5.9]; % remainder of X-,along-canyon survey
    case 'wawaves14'
        DelDNtoBT = -0.44;
        DelDNtoUP = 0.47; % SW = 09000
        MagDcl = 17; % JM matlab magdev, will use on TRBM also
        EnsDNtoUP = [238.72,229; 240.87,0; 241.4798,144];
        DelDNtoSW = ...
            [238.73,2.5; 238.76,2.2; % test off Oregon coast
            240.87,3.0; 241.34,1.0; % first wave chasing
            241.4798,0.8; 241.88,0.0; % back in after mmp
            241.8801,4.0; 241.986,2.5; % changed clock during short transit, then
            242.0265,4.0; 242.55, 2.5; % after changing SU uir-file, deeper casts
            243.0,1.3; 243.5,-0.3; 243.917,-1.7]; % end of cruise
    case 'arcticmix'
        DelDNtoBT = -0.44;
        DelDNtoUP = 0.47; % SW = 09000
        MagDcl = 21.0; % geomag, and magdev.m @ 72.49Nx144.97W
        EnsDNtoUP = [...
            238.9530,0;...
            243.2510,35;...
            244.9960,144;...
            245.821,0;...
            247.114,2;...
            249.8950,-40;...
            253.737,283;...
            255.828,34;
            256.819,0;...
            257.301,39;...
            260.459,0];
        DelDNtoSW = ...
            [238.9530,56.538;...
            243.251,3.951;...
            243.5,4.32;...
            243.75,4.72;...
            243.903,4.954;...
            244.996,7.163;...
            245.660,8.186;...
            245.821,8.478;...
            247.114,10.596;...
            249.895,2.714;...
            250.576,3.721;...
            253.737,23.602;...
            254.398,24.607;...
            255.054,25.617;...
            255.7,26.618;...
            255.828,1.482;...
            256.819,3.179;...
            257.301,3.915;...
            257.503,4.186;...
            258.196,5.192;...
            260.459,10.165;...
            261.14,11.165]; % sample data for dave
    case 'sproultest'
        DelDNtoBT = -0.44;
        DelDNtoUP = 0.47; % SW = 09000
        MagDcl = 11.73;
        EnsDNtoUP = [76.5 36; 
                    76.81 0;
                    77.07 22;
                    77.16 40;
                    77.48 27;
                    77.65 37;
                    77.87 0;
                    78.5  29756
                    78.51 30336;
                    78.59 325;
                    78.6  0]; % DtoU
        DelDNtoSW = [76.5 -0.2860; 
                    76.78 -0.2430;
                    76.81 0.02;
                    76.82 -0.235;
                    76.85 -0.2040;
                    76.9 -0.174;
                    77.07 -0.0810;
                    77.16 0.008;
                    77.3 0.094;
                    77.48 0.211;
                    77.65 0.359;
                    77.8 0.4310;
                    77.87 0.499;
                    78 0.533;
                    78.2 0.65;
                    78.5 0.83;
                    78.6 0.921]; % DNtoCTD
        DelUPtoSW = [76.5 1.8370; 
                    76.78 1.8280;
                    76.81 1.9350;
                    76.82 1.7540;
                    76.85 1.717;
                    76.9  1.672;
                    77.07 1.45;
                    77.16 1.344;
                    77.3 1.193;
                    77.48 0.949;
                    77.65 0.77;
                    77.8 0.566;
                    77.87 0.47;
                    78 0.311;
                    78.2 0.065;
                    78.5 -0.318;
                    78.6 -0.449]; % UPtoCTD
    case 'fleat'
        DelDNtoBT = -0.44;
        DelDNtoUP = 0.47;
        MagDcl = 17/60; % geomag, and magdev.m @ 72.49Nx144.97W
        EnsDNtoUP = [...
            156,29;...
            160.656,-1;...
            161.65,-2];
        DelDNtoSW = ...
            [156,4.0890;...
            156.9956, 3.7;...
            157.01,3.6;...
            157.03,3.5;...
            157.05,3.4;...
            157.09,3.2;...
            157.12,3.0;...
            157.16,2.8;...
            157.19,2.6;...
            157.23,2.4;...
            157.27,2.2;...
            157.31,2.0;...
            157.34,1.8;...
            157.38,1.6;...
            157.42,1.4;...
            157.455,6.021;...
            157.49,5.8;...
            157.53,5.6;...
            157.56,5.4;...
            157.6,5.2;...
            157.64,5;...
            157.67,4.8;...
            157.74,4.4;...
            157.81,4;...
            157.832,5.421;...
            157.86,5.283;...
            157.9,5.0690;...
            157.98,4.635;...
            158.05,4.2560;...
            158.15,3.7180;...
            158.2,3.43;...
            160.656,3.746;...
            160.7,3.438;...
            160.8,2.94;...
            160.9,2.42;...
            161.1,1.408;...
            161.158,3.789;...
            161.2,3.737;...
            161.65,3.118;...
            161.7,3.083;...
            161.8,2.676;...
            161.9,2.446;...
            162,2.352;...
            162.1,2.421;...
            162.2,2.265;...
            ];
    case 'lajit2'
        DelDNtoBT = -0.44;
        DelDNtoUP = 0.47; % SW = 09000
        MagDcl = 11.73;
        EnsDNtoUP = [...
            253.68, -2;...
            253.772,-3;...
            254.24,-5];
        DelDNtoSW = ...
            [253.68,3.8360;...
            253.70,3.8240;...
            253.753,3.557;...
            253.755,3.561;...
            253.772,3.32;...
            253.774,3.302;...
            254.24,2.740;...
            254.4,2.642;...
            254.6,2.453;...
            254.8,2.249;...
            255, 2.0640;...
            255.5, 1.576;...
            ];
    otherwise
        DelDNtoBT = -0.44;
        DelDNtoUP = 0.47; % SW = 09000
        MagDcl = 19.5 ; % for Puget Sound, Hood Canal, approx
        % for no BTing after 2003:
        % DelDNtoBT = -0.1; % doesn't matter, not used
        % DelDNtoUP = 0.61;
end

% convert to yearday increments
DelDNtoBT = DelDNtoBT / 86400;
DelDNtoUP = DelDNtoUP / 86400;
DelDNtoSW(:,2) = DelDNtoSW(:,2) / 86400;
if ~isempty(DelUPtoSW)
    DelUPtoSW(:,2) = DelUPtoSW(:,2) / 86400;
end