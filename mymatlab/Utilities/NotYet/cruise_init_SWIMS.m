% cruise_init_SWIMS.m  -  initialize folders, indices for cruise;
%	Edit values above line='%%%%%%%%%%%%%%%%%% HERE %%%%%%%%%%%%%%%%%%'
%   for new cruise.  cd to X:/swims/'cruisename' before running this.
% Dave W, apr-2002

CR_name = 'ps02';
SWIM_no = 2;
% structures for matlab index files:
Cruise.name = 'ps02';
Cruise.start_date = '25-Apr-2002 10:00:00';  % UTC, approx to start
Cruise.end_date = '11-May-2002 08:00:00';
Def(1).year = 2002;
Def(1).start_yday = 0;
Def(1).end_yday = 400;
Def(1).yday_seconds_offset = 0;
Idx(1).yday_beg = [];
Idx(1).yday_end = [];
Idx(1).filename = [];
PROG ='cruise_init_SWIMS.m';

cr_dir=pwd;
DATA_dir = 'C:/swims/data/';

%%%%%%%%%%%%%%%%%% HERE %%%%%%%%%%%%%%%%%%

data_typs={'CTD';'GPS';'LD';'TD';'ADCP';'ADUP';'ADDN';'KVH'};
cr_flds = {'indexes';'griddata';'data_mat'};

disp(['SET UP for new cruise=' CR_name ' in: ' cr_dir '  -  ']);
disp(['(1) Folders = ' cr_flds'])
disp(['  and subfolders in ./data_mat and ' DATA_dir ' for: ')
disp(data_typs'), disp('')

x = input('Ready to proceed? ', 's');
switch x
case {'yes'}
   disp('Okay, starting ...');
   more = 1;
otherwise
   disp('Skipping folder creation !');
   more = 0;
end

%% Create new folders
if more
    for id=1:length(cr_flds)
        [sta,msg] = mkdir(cr_dir, cr_flds{id});
        if ~sta, disp(msg); end
    end
    for id=1:length(data_typs)
        [sta,msg] = mkdir(fullfile(cr_dir,'data_mat'), data_typs{id});
        if ~sta, disp(msg); end
        [sta,msg] = mkdir(DATA_dir, data_typs{id}); % will often already exist
    end
    disp('  Done making folders.')   
end

disp(''),disp('(2) Initialize indices for managing raw data files;')
disp(['    create in ' fullfile(cr_dir,'indexes')]), disp('')

x = input('Continue? ', 's');
switch x
case {'yes'}
   disp('Okay, starting ...');
   more = 1;
otherwise
   disp('Skipping index creation - THESE WILL BE NEEDED !');
   more = 0;
end

%% initialize raw,matlab index files
if more
    % raw index first
    file_name=[]; file_yday=[]; file_status=[]; file_copied=[];
    for id=1:length(data_typs)
        fnidx = [data_typs{id} '_' CR_name '_rawfiles.mat']
        disp(['Init ' fnidx])
        if exist(fullfile(cr_dir,'indexes',fnidx), 'file')
            disp('   Already exists, skipping !!')
        else
            save( fullfile(cr_dir,'indexes',fnidx), ...
                'file_name','file_yday','file_status','file_copied');
        end
    end % of initializing index into raw (acquisition) files
    
    % then matlab (stage 1) index; 
    for id=1:length(data_typs)
        fnidx = [data_typs{id} '_' CR_name '_matfiles.mat']
        disp(['Init ' fnidx])
        if exist(fullfile(cr_dir,'indexes',fnidx), 'file')
            disp('   Already exists, skipping !!')
        else
            clear Set_params Index
            Set_params = Def;  Index = Idx;
            % add more fields to Set_params(1), Index(1) if needed
            %
            save( fullfile(cr_dir,'indexes',fnidx), ...
                'Cruise','Set_params','Index','PROG');
        end
    end % of initializing index into Matlab (stage 1) files

end
