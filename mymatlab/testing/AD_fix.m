%Fix ADDN, ADUP mat-files with out-of-order entries (index okay, by accident)

DD='D:\swims\BS03\data_mat\ADDN';
FF = 'ADDN-2003-074-1705-01.mat';

load([DD '\' FF])

iok = find(ADDN.yday > 74);

fnms = fieldnames(ADDN);

for i=1:length(fnms)
    switch fnms{i}
    case 'z_adcp'
        x=1;
    otherwise
        eval(['ADDN.' fnms{i} '= ADDN.' fnms{i} '(:,iok);']);
    end
end

DD='D:\swims\BS03\data_mat\ADUP';
FF = 'ADUP-2003-074-1705-01.mat';

load([DD '\' FF])

iok = find(ADUP.yday > 74);

fnms = fieldnames(ADUP);

for i=1:length(fnms)
    switch fnms{i}
    case 'z_adcp'
        x=1;
    otherwise
        eval(['ADUP.' fnms{i} '= ADUP.' fnms{i} '(:,iok);']);
    end
end
