% Run_ADCP_EnsYdays.m - parse SWIMS ADCP ensemble timestamp files, and save
%  in one large matlab data file for each type (ADCP1 and ADCP2) dpw - 9/15

cruise='SproulTest';
RawFld = fullfile(matlabdisk,'swims','SproulTest','raw_data');
%RawFld = fullfile('C:','swims','data'); % ADCP*_ensembles_logfile.txt are here
MatFld = swimsmatdata; % file both here

RawTyps = {'ADCP1','ADCP2'};
MatTyps = {'ADDN','ADUP'};

for ityp = 1:2
    MatFil = fullfile(MatFld, [MatTyps{ityp} '-EnsYdays.mat']);
    RawFils = dir( ...
        fullfile(RawFld, [RawTyps{ityp} '_*_ensembles_logfile.txt']) );
    % assuming filenames follow chronological order
    
    nFilB = 1; % by default, start with file # 1
    LastPCyday = -1;
    if exist(MatFil, 'file');
        load(MatFil)
    else % initialize
        EnsYdays.ens_no = [];
        EnsYdays.RDyday = [];
        EnsYdays.PCyday = [];
        EnsYdays.LastRaw = [];
    end
    
    % re-do last previous file, and those that follow
    if ~isempty(EnsYdays.LastRaw)
        LastPCyday = EnsYdays.PCyday(end);
        for i=1:length(RawFils)
            if strcmp(RawFils(i).name, EnsYdays.LastRaw)
                nFilB = i;
                break; % start with this file
            end
        end
    end
    
    % Loop thru specified raw files
    for ifil = nFilB:length(RawFils)
        if RawFils(ifil).bytes > 20
            disp([' processing ' RawFils(ifil).name ' --'])
            EY = Read_rawADCP_EnsYday(fullfile(RawFld, RawFils(ifil).name));
            if isempty(EY)
                disp('    -- no valid data')
                continue; % no valid data in this file
            end
            % skip portion (of first file) that was already processed
            iok = find(EY.PCyday > LastPCyday);
            if ~isempty(iok) % append to existing data
                EnsYdays.ens_no = [EnsYdays.ens_no EY.ens_no(iok)];
                EnsYdays.RDyday = [EnsYdays.RDyday EY.RDyday(iok)];
                EnsYdays.PCyday = [EnsYdays.PCyday EY.PCyday(iok)];
                EnsYdays.LastRaw = RawFils(ifil).name;
            else
                disp('    -- no new data')
            end
        end
    end
    
    % save updated matlab structure
    if ~isempty(EnsYdays.PCyday)
        disp(['Saving updated EnsYdays.* in ' MatFil])
        save(MatFil, 'EnsYdays')
    end
end