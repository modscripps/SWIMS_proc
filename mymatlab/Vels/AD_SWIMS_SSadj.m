% AD_SWIMS_SSadj.m - with all set up, compute velocities: result in ADX
%  This version will use axes and transforms specified by RDI -- unlike MHA's and DPW's, does
%  account for pitch vs. roll rotations of axes (i.e. w is supposed to be straight up, u,v horiz)
%  This version (twist) maps UPward beams and pitch/roll onto DNward, to use same txforms.
%% Version V2 tries variations on error checking:

%% Place results in structure ADX (initialize)
% VelDN.*, VelUP.*, en_b,en_e prepared in AD_SWIMS_PrSynch.m 
clear ADX
%ADX.depth = [6:2:600]'; Trow = length(ADX.depth); % Mamala Bay
%ADX.depth = [6:2:660]'; Trow = length(ADX.depth); % Kauai Ridge
% ADX.depth = [6:2:250]'; Trow = length(ADX.depth); % ML04
% ADX.depth = [6:2:250]'; Trow = length(ADX.depth); % stf07
% ADX.depth = [6:2:600]'; Trow = length(ADX.depth); % philex08, later
% ADX.depth = [6:2:250]'; Trow = length(ADX.depth); % hc03/ps02
% ADX.depth = [6:2:150]'; Trow = length(ADX.depth); % mc09
% ADX.depth = [6:2:150]'; Trow = length(ADX.depth); % MORT
% ADX.depth = [6:2:580]'; Trow = length(ADX.depth); % wa_nliw_apr2013
ADX.depth = [4:2:580]'; Trow = length(ADX.depth); % WaWaves14, ArcticMix
% Map pings/ensembles to columns - depends on synch/asynch
icDN = iDN*NaN;
icUP = iUP*NaN;
if SynchUPDN || length(iUP) < 2
    % (Use Synch--via ~isnan(EnsDNtoUP)--even if no UP data are available)
    e1 = min(en_bD, en_bU);
    e2 = max(en_eD, en_eU);
    ncol = e2-e1 + 1;
    icDN = VelDN.ens_no(iDN) - e1 + 1;
    icUP = VelUP.ens_no(iUP) - e1 + 1;
elseif length(iDN) < 2
    % UP data only, possibly asynch: use ensemble numbers
    ncol = length(iUP);
    icUP = 1:ncol;
else % Asynch, data from DN and UP
    % For duration of type (UP/DN) with shorter intervals, use that one's
    % pings with closest pings from other type (which skips some columns).
    % Pre/Append pings from slower type outside faster one's range.
    ybD = VelDN.yday_DN(iDN(1));
    yeD = VelDN.yday_DN(iDN(end));
    dtD = median( diff(VelDN.yday_DN(iDN)) );
    ydD = [ybD:dtD:yeD];
    ybU = VelUP.yday_UP(iUP(1));
    yeU = VelUP.yday_UP(iUP(end));
    dtU = median( diff(VelUP.yday_UP(iUP)) );
    ydU = [ybU:dtU:yeU];
    % pick finer yday grid, if indicated add coarser one on ends
    if dtU>dtD
        ydM = ydD; ydE = ydU;
    else
        ydM = ydU; ydE = ydD;
    end
    ix = find(ydE<ydM(1));
    ydM = [ydE(ix) ydM];
    ix = find(ydE>ydM(end));
    ydM = [ydM ydE(ix)];
    % map ensembles
    ncol = length(ydM);
    icDN = interp1(ydM, [1:ncol], VelDN.yday_DN(iDN), 'nearest');
    icUP = interp1(ydM, [1:ncol], VelUP.yday_UP(iUP), 'nearest');
    % any duplicate column numbers will be skipped during processing
    clear ybD yeD ydD ybU yeU ydU ydM ydE
end

% initialize arrays with NaN's
ADX.u = NaN * ones(length(ADX.depth), ncol);
ADX.v = ADX.u; ADX.w = ADX.u; ADX.werr = ADX.u;
% More, for diag data screening
ADX.ec1=ADX.u;ADX.ec2=ADX.u;ADX.ec3=ADX.u;ADX.ec4=ADX.u;
% initialize vectors
Avars = {'yday_DN','yday_UP','yday_BT','pr_DN','pr_UP','pr_BT','bottomBT','botBTmin', ...
    'ens_DN','ens_UP','headT_DN','headT_UP','headT','uBT','vBT','wBT','werrBT', ...
    'roll_DN','pitch_DN','roll_UP','pitch_UP','w_DN','w_UP'};
for iv=1:length(Avars)
    eval(['ADX.' Avars{iv} '=ADX.u(1,:);']);
end

%Now set up constants and parameters
SvelA = 1500; % Used in ADCPs' setups
dz_DN = 0.21; dz_UP = -0.21; % offsets to txducers (17" total)

% Fill in depth vs time arrays, VelDN first    
for ip=1:length(iDN)
    wh = iDN(ip); ieN = VelDN.ens_no(wh); icol = icDN(ip);
    if isnan(icol)
        continue; % ensemble didn't map to any column
    end
    if ~isnan(ADX.ens_DN(icol))
        continue; % duplicate for column, use only first one
    end
    ADX.ens_DN(icol) = ieN;

    % The first step is to compute depth vectors adjusted for sound speed
    % (approx, using value at CTD only).
    SSadj = VelDN.svel_DN(wh)/SvelA;
    z_rel = VelDN.z_adcp; z_rel(1) = VelDN.z_adcp(1) * SSadj;
    z_rel(2:end) = z_rel(1) + cumsum([diff(VelDN.z_adcp)*SSadj]);
    %Now the absolute depth vector is obtained by adding on SWIMS's depth at the time of the ping.
    %z_dep = VelDN.pr_DN(wh)*100 + z_rel + dz_DN;
    z_dep = VelDN.pr_DN(wh) + z_rel + dz_DN; % pr/dbar
    %Now correct the bottom ranges for sound speed  
    SSadjBT = VelDN.svel_BT(wh)/SvelA;
    rangeBT = VelDN.rangeBT(wh) * SSadjBT;
    rangeBTmin = VelDN.rangeBTmin(wh) * SSadjBT;
    
    % Now interpolate velocity components and echo amplitudes onto the grid.
    % Make sound speed adjustment to velocities here, also.
    vgs = {'u','v','w','werr', 'ec1','ec2','ec3','ec4'};
    for i=1:8
        vv = VelDN.(vgs{i})(:,wh);
        if i<5
            vv = vv*SSadj; % for velocity components
        end
        ig = find(~isnan(vv));
        if length(ig)>1
            gv = interp1(z_dep(ig(1):ig(end)), vv(ig(1):ig(end)), ADX.depth);
            ik = find(~isnan(gv)); % grid depths of good data
            ADX.(vgs{i})(ik,icol) = gv(ik); % insert into output grid
        end
    end
    
    % Now update vectors in the output structure.
    %ADX.bottomBT(icol) = rangeBT + VelDN.pr_BT(wh)*100 + dz_DN; % bottom depth
    ADX.bottomBT(icol) = rangeBT + VelDN.pr_BT(wh) + dz_DN; % pr/dbar
    ADX.botBTmin(icol) = rangeBTmin; % minimum vertical range, VelDN-to-bottom
    % sound speed adjustments to BT velocity
    ADX.uBT(icol) = VelDN.uBT(wh)* SSadjBT;
    ADX.vBT(icol) = VelDN.vBT(wh)* SSadjBT;
    ADX.wBT(icol) = VelDN.wBT(wh)* SSadjBT;
    ADX.werrBT(icol) = VelDN.werrBT(wh)* SSadjBT;

    Avars = {'yday_DN','yday_BT','pr_DN','pr_BT','w_DN'};
    for iv=1:length(Avars)
        ADX.(Avars{iv})(icol) = VelDN.(Avars{iv})(wh);
    end
    ADX.pitch_DN(icol) = VelDN.pitch(wh);
    ADX.roll_DN(icol) = VelDN.roll(wh);
    ADX.headT_DN(icol) = VelDN.headT(wh);
end  % of VelDN ensemble processing (with BT's)

% Fill in depth vs time arrays, VelUP second   
for ip=1:length(iUP)
    wh = iUP(ip); ieN = VelUP.ens_no(wh); icol = icUP(ip);
    if isnan(icol)
        continue; % ensemble didn't map to any column
    end
    if ~isnan(ADX.ens_UP(icol))
        continue; % duplicate for column, use only first one
    end
    ADX.ens_UP(icol) = ieN;
    
    % The first step is to compute depth vectors adjusted for sound speed
    % (approx, using value at CTD only).
    SSadj = VelUP.svel_UP(wh)/SvelA;
    z_rel = VelUP.z_adcp; z_rel(1) = VelUP.z_adcp(1) * SSadj;
    z_rel(2:end) = z_rel(1) + cumsum([diff(VelUP.z_adcp)*SSadj]);
    %Now the absolute depth vector is obtained by adding on SWIMS's depth at the time of the ping.
    %z_dep = VelUP.pr_UP(wh)*100 + z_rel + dz_UP;
    z_dep = VelUP.pr_UP(wh) + z_rel + dz_UP; % pr/dbar
    
    % Now interpolate velocity components and echo amplitudes onto the grid.
    % Make sound speed adjustment to velocities here, also.
    vgs = {'u','v','w','werr', 'ec1','ec2','ec3','ec4'};
    for i=1:8
        vv = VelUP.(vgs{i})(:,wh);
        if i<5
            vv = vv*SSadj; % for velocity components
        end
        ig = find(~isnan(vv));
        if length(ig)>1
            gv = interp1(z_dep(ig(1):ig(end)), vv(ig(1):ig(end)), ADX.depth);
            ik = find(~isnan(gv)); % grid depths of good data
            ADX.(vgs{i})(ik,icol) = gv(ik); % insert into output grid
        end
    end
    
    % Now update vectors in the output structure.
    Avars = {'yday_UP','pr_UP','w_UP'};
    for iv=1:length(Avars)
        ADX.(Avars{iv})(icol) = VelUP.(Avars{iv})(wh);
    end
    ADX.pitch_UP(icol) = VelUP.pitch(wh);
    ADX.roll_UP(icol) = VelUP.roll(wh);
    ADX.headT_UP(icol) = VelUP.headT(wh);
end % of VelUP ensemble processing

% return
ADX_check % to apply screening algorithms

%keyboard

% rotate to SWIMS coordinates: facing nose (front-to-back):
% along(v, pos flow fwd, usually < 0), across(u, pos flow to left) SWIMS
c45=cos(45*pi/180); s45=sin(45*pi/180);
% VEL_rel = uDN x_hat' + vDN y_hat' = uSW x_hat + vSW y_hat
% So both (VelDN,SWIMS) are RH axes from above; angle from x_hat' to x_hat = 45
clear xx yy
xx = ADX.u; yy = ADX.v;
ADX.u = xx*c45 + yy*s45;
ADX.v = -xx*s45 + yy*c45;
xx = ADX.uBT; yy = ADX.vBT;
ADX.uBT = xx*c45 + yy*s45;
ADX.vBT = -xx*s45 + yy*c45;

% Use average of DN,UP true heading if both are non-NaN;
%  watch for one slightly NNE>0, other NNW<360
%******Remember heading is not corrected at present for roll/pitch effects.
ibu = find(isnan(ADX.headT_UP));
ibd = find(isnan(ADX.headT_DN));
ibN = find((ADX.headT_UP>300&ADX.headT_DN<60) | (ADX.headT_UP<60&ADX.headT_DN>300));
ADX.headT = (ADX.headT_UP+ADX.headT_DN)/2;
ADX.headT(ibu) = ADX.headT_DN(ibu); 
ADX.headT(ibd) = ADX.headT_UP(ibd); 
ADX.headT(ibN) = ADX.headT(ibN)+180; 

%keyboard

Ang = -(360-ADX.headT)*pi/180; % math convention (here, SWIMS to E,N)
for i=1:length(Ang)
    xx = ADX.u(:,i); yy = ADX.v(:,i);
    ADX.u(:,i) = xx*cos(Ang(i)) + yy*sin(Ang(i));
    ADX.v(:,i) = -xx*sin(Ang(i)) + yy*cos(Ang(i));
    xx = ADX.uBT(i); yy = ADX.vBT(i);
    ADX.uBT(i) = xx*cos(Ang(i)) + yy*sin(Ang(i));
    ADX.vBT(i) = -xx*sin(Ang(i)) + yy*cos(Ang(i));
end

