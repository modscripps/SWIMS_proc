function LinLims = find_profiles_LineOut(LDat);

LinLims = []; 
LnDel = 5 * 10; % number of LineOut points to compute over [5s*10Hz]
%LnDel = 8 * 0.5; % number of LineOut points to compute over [0.5-Hz data]
% Make sure enough data are supplied
if length(LDat.line_out) < LnDel*5
    return
end

%Ldt = LDat.eltim_out(LnDel+1:end) - LDat.eltim_out(1:end-LnDel);
Ldt = (LDat.PCyday(LnDel+1:end) - LDat.PCyday(1:end-LnDel)) * 86400;
Ldl = LDat.line_out(LnDel+1:end) - LDat.line_out(1:end-LnDel);
Lrat = Ldl./(Ldt/60); % in m/min
Ldy = (LDat.PCyday(LnDel+1:end) + LDat.PCyday(1:end-LnDel)) / 2;
ix = find(Ldt<=0 | Ldt>10*LnDel); % probably not continuous data; gap/delay
Lrat(ix) = 0;
Ldir = sign(Lrat); % neg=in/up; pos=out/down; zero=pause (switches sense later)
Lirev = find( abs(diff(Ldir)) >0 );
%keyboard
LCydb = []; LCyde = []; LCdir = [];
for ir=1:length(Lirev)
    ip = Lirev(ir); ic = ip+1; % previous and current index
    if Ldir(ip)~=0 % finish previous up/down
        LCyde(end+1) = Ldy(ip);
        if ir==1 % start and type of first one
            LCydb(1) = Ldy(1);
            LCdir(1) = Ldir(ip);
        end
    end
    if Ldir(ic)~=0 % start and type of current one
        LCydb(end+1) = Ldy(ic);
        LCdir(end+1) = Ldir(ic);
        if ir==length(Lirev) % finish last one
            LCyde(end+1) = Ldy(end);
        end
    end
end
ix = find( (LCyde-LCydb)*86400 < 5); % Too short in duration
LCydb(ix)=[]; LCyde(ix)=[]; LCdir(ix)=[];
% clear LDat Ldt Ldl Ldy Ldir Lirev Lrat
% % assume 1.5-s delay in Line data timestamps (APL, relative to CTD scans)
% LinLims.yday_beg = LCydb - 1.5/86400;
% LinLims.yday_end = LCyde - 1.5/86400;

% NO delay in Line data timestamps from LCI90 (relative to CTD scans)
LinLims.yday_beg = LCydb;
LinLims.yday_end = LCyde;

LinLims.type = -LCdir; % switch to up>0, down<0
