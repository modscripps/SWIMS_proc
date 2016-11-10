%Fix ADDN, ADUP mat-files with out-of-order entries (index okay, by accident)

ADtyp={'DN','UP'};

clear ADinfo
for id=1:2
    DD = ['D:\swims\BS03\data_mat\AD' ADtyp{id} '\']
    Mfs = dir([DD 'AD' ADtyp{id} '-*']);
    IDx = ['D:\swims\BS03\indexes\AD' ADtyp{id} '_BS03_matfiles.mat']
    % Cruise      Index       PROG        Set_params
    % Index       yday_beg yday_end filename: {1x28 cell}
    load(IDx)
    ydb = []; yde = []; ydL = 69; clear fn
    ADinfo(id).ens=[]; ADinfo(id).yday=[]; ADinfo(id).fno=[];ADinfo(id).bot=[];
    figure(id)
    for ix=1:length(Mfs)
        fn{ix} = Mfs(ix).name;
        clear ADDN ADUP
        load([DD fn{ix}])
        eval(['ydy = AD' ADtyp{id} '.yday;']);
        ydb(ix) = min(ydy); yde(ix) = max(ydy);
        if ~isempty(find(diff([ydL ydy])<=0))
            disp([fn{ix} ' not monotonic'])
        else
            disp(fn{ix})
        end
        eval(['ADinfo(id).ens=[ADinfo(id).ens AD' ADtyp{id} '.ens_no];']);
        ADinfo(id).yday=[ADinfo(id).yday ydy];
        ADinfo(id).fno=[ADinfo(id).fno ix*ones(size(ydy))];
        if id<2
            eval(['ADinfo(id).bot=[ADinfo(id).bot AD' ADtyp{id} '. bottomBT];']);
        end
        c='rg'; c=c(mod(ix,2)+1);
        plot(ydy, diff([ydL ydy])*86400, [c '-'], ...
            ydy, diff([ydL ydy])*86400, [c '.']), hold on
        ydL = ydy(end);
    end
    clear Index
    Index.yday_beg = ydb; Index.yday_end = yde;
    Index.filename = fn;
    save(['./AD' ADtyp{id} '_BS03_matfiles.mat'], 'Cruise', 'Index', 'PROG', 'Set_params')
    keyboard
end