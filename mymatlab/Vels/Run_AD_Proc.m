function out = Run_AD_Proc()

% Run_AD_Proc.m  -  loop thru ADDN,ADUP indexes, create data files with
% velocity components u,v,w,werr (and for BT, if present) such that:
%% positive u is flow from NOSE RIGHT to TAIL LEFT (from DN beam 1 to 2),
%% positive v is flow from TAIL RIGHT to NOSE LEFT (from DN beam 4 to 3),
%% positive w is flow towards ADDN (SWIMS descending).
% These are in the leveled frame (pitch, roll taken out), but depths/ranges
% and velocity amplitudes are not yet adjusted for sound speed variations
% from the nominal 1500 m/s (This will be done in subsequent processing,
% as data are mapped to actual depths by synching with SWIMS CTD pressure)

cruise='SproulTest'; 
TTYYPP = {'UP'; 'DN'};
TopFld = '/Users/ecfine/Documents/MATLAB/swims/SproulTest';

for iittyypp = 1:2
    ADtyp = ['AD' TTYYPP{iittyypp}];
    VLtyp = ['Vel' TTYYPP{iittyypp}];
    ADidx = fullfile(TopFld,'indexes', [ADtyp '_' cruise '_matfiles.mat']);
    ADfld = fullfile(TopFld,'data_mat', ADtyp);
    ADSTR = load(ADidx);
    clear Cruise Set_params PROG Index
    
    % Initialize new index for pre-processed velocity
    VLidx = fullfile(TopFld,'indexes', [VLtyp '_' cruise '_matfiles.mat']);
    VLfld = fullfile(TopFld,'data_mat', VLtyp);
    
    iSet = 1;
    % Undo adjustments to recorded yeardays, if they have been defined
    % they are corrected here (WATCH OUT: make sure this is done correctly!!)
    SpS = ADSTR.Set_params(iSet);
    YD_jump_aft = [];
    YD_jump_by = [];
    if isfield(SpS, 'yday_jump_after')
        YD_jump_aft = [SpS.yday_jump_after,  inf]; % last one to adjust to end
        YD_jump_by = SpS.yday_jump_adjust;
    end
    YD_off = SpS.yday_seconds_offset / 86400;

    % In case of pre-existing adjustments in VelUP/DN, use them here
    YV_jump_aft = [];
    YV_jump_by = [];
    YV_off = 0;
    if exist(VLidx, 'file')
        VLSTR = load(VLidx);
        if isfield(VLSTR.Set_params, 'yday_jump_after')
            YV_jump_aft = ...
                [VLSTR.Set_params.yday_jump_after,  inf]; % last one to adjust to end
            YV_jump_by = VLSTR.Set_params.yday_jump_adjust;
        end
        YV_off = VLSTR.Set_params.yday_seconds_offset / 86400;
        % Now, find yearday of last pre-processed ADDN/UP ensemble
        yb = VLSTR.Index(iSet).yday_beg;
        ye = VLSTR.Index(iSet).yday_end;
        % Adjust limits in manner of get_ADupdn_data.m
        for i=1:length(YV_jump_by)
            ix = find( yb>=YV_jump_aft(i) & yb<YV_jump_aft(i+1) );
            yb(ix) = yb(ix) + YV_jump_by(i);
            ix = find( ye>=YV_jump_aft(i) & ye<YV_jump_aft(i+1) );
            ye(ix) = ye(ix) + YV_jump_by(i);
        end
        yd_b = yb + YV_off; yd_e = ye + YV_off;
        YDnew_beg = yd_e(end);
    else % initialize new index
        clear VLSTR
        VLSTR.Cruise = ADSTR.Cruise(iSet);
        VLSTR.Set_params = ADSTR.Set_params(iSet);
        VLSTR.PROG = 'Run_AD_Proc.m';
        VLSTR.Index.yday_beg = []; VLSTR.Index.yday_end = [];
        VLSTR.Index.filename = [];
        if isfield(SpS, 'yday_jump_after')
            VLSTR.Set_params.yday_jump_adjust = 0 * SpS.yday_jump_adjust;
        end
        VLSTR.Set_params.yday_seconds_offset = 0;
        YDnew_beg = 0;
    end
    
    %% Now, check begin,end yeardays (adjusted) of ADDN/UP matlab files
    %%  to decide which ones need to be processed into VelDN/UP files.
    %% Start after last one completed, and stop with next to the last
    %%  one (in case the last one is still being added to)
    yb = ADSTR.Index(iSet).yday_beg;
    ye = ADSTR.Index(iSet).yday_end;
    % Adjust limits in manner of get_ADupdn_data.m
    for i=1:length(YD_jump_by)
        ix = find( yb>=YD_jump_aft(i) & yb<YD_jump_aft(i+1) );
        yb(ix) = yb(ix) + YD_jump_by(i);
        ix = find( ye>=YD_jump_aft(i) & ye<YD_jump_aft(i+1) );
        ye(ix) = ye(ix) + YD_jump_by(i);
    end
    yd_b = yb + YD_off; yd_e = ye + YD_off;
    iFb = find(yd_b > YDnew_beg); 
    if isempty(iFb)
        continue;
    end
    iFb = iFb(1); % next ADDN/UP after last VelDN/UP
    
    iFe = length(yd_b)-1; % next to last ADDN/UP matlab file
    
    %%% Uncomment the next line at end of cruise, when all AD*.mat files are complete
    iFe = iFe+1; disp('NOTE: doing last ADCP file!!') 
    
    if iFb>iFe
        disp(['No ' VLtyp ' files ready for processing.'])
        continue
    end
    disp(['Create ' VLtyp ' files (and index entries) for:'])
    disp([ADSTR.Index(iSet).filename{iFb} ' thru ' ADSTR.Index(iSet).filename{iFe}])
    x = input('Okay to proceed (existing ones safe somewhere) ? Type "Y" if so!', 's');
    if isempty(x) || ~strcmp(x(1),'Y')
        disp('NOT STARTED!')
        continue
    end
    disp('Okay - starting ...')
    if ~isempty(YV_jump_by) % New VelDN/UP will have correct time, turn off 'jumps'
        if abs(YV_jump_by(end)) > 1e-5
            VLSTR.Set_params.yday_jump_adjust(end+1) = 0;
            VLSTR.Set_params.yday_jump_after(end+1) = YDnew_beg;
        end
    end

    % Loop thru ADDN/ADUP files, generate and save preprocessed file
    for iFil = iFb:iFe
        Fnam = ADSTR.Index(iSet).filename{iFil};
        VLFnam = [VLtyp Fnam(5:end)];
        yb = ADSTR.Index(iSet).yday_beg(iFil);
        ye = ADSTR.Index(iSet).yday_end(iFil);
        % Adjust limits in manner of get_ADupdn_data.m
        for i=1:length(YD_jump_by)
            ix = find( yb>=YD_jump_aft(i) & yb<YD_jump_aft(i+1) );
            yb(ix) = yb(ix) + YD_jump_by(i);
            ix = find( ye>=YD_jump_aft(i) & ye<YD_jump_aft(i+1) );
            ye(ix) = ye(ix) + YD_jump_by(i);
        end
        yd_b = yb + YD_off; yd_e = ye + YD_off;
        %
        if strcmp(cruise,'bs03') && strcmp(ADtyp,'ADUP') && yd_b < 87.0
            continue; % ADUP velocity data were worthless before this
        end
        
        % Compute pre-processed velocity, etc; Save data file, update index
        disp(['Processing ' VLFnam])
        AD_SWIMS_Process
        clear(VLtyp)
        eval([VLtyp ' = ADP;'])
        save( fullfile(VLfld, VLFnam), VLtyp )
        VLSTR.Index.filename{end+1} = VLFnam;
        VLSTR.Index.yday_beg(end+1) = min(ADP.yday);
        VLSTR.Index.yday_end(end+1) = max(ADP.yday);
        Index = VLSTR.Index; Set_params = VLSTR.Set_params;
        Cruise = VLSTR.Cruise; PROG = VLSTR.PROG;
        save(VLidx, 'Set_params', 'Index', 'Cruise', 'PROG')
        pause(0.2)
    end
        
    disp('DONE!')
end
end
