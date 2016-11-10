function [AllOverturns,CTDout]=FindOverturnsFCN(CTD,PP);
%function [AllOverturns,CTDout]=FindOverturnsFCN(CTD,PP);
%Do thorpe sorting on a CTD profile or array of CTD profiles.  Pass in CTD, a structure with the
%fields z, T, S, H, lat, lon, yday and an id.  The names of the fields may
%be specified in the params structure PP.
%
%The calculation occurs in two steps, of which this file is the first.
%
%First, the routine performs basic
%checking for NaN's and pressure reversals and then calls the function
%OverturnListFCN which returns a list of Overturn Structures.  The routine uses the cumulative Thorpe
%displacement to find all reordering regions.  Each of
%these structures has the basic quantities such as rms Thorpe displacement,
%start and end depths, etc, but also contains a structure, err, obtained by
%calling the function GKerr (see help file for that).  This structure
%contains the expected turbulent and noise signals in the stratification,
%N, of the overturn.  Noise signals result from both depth and density
%precision are computed.  Also, the expected Thorpe displacement for K_rho
%= 10^-3 and 10^-4 are computed.  This structure allows observed quantities
%to be compared with those expected from noise and turbulence.
%
%Second, the function EpsFromOverturns2 is called, which takes the
%overturns, checks them according to the tests and limits specified in PP,
%and computes mean dissipation rate by storing that for each overturn in
%the depth range occupied by that overturn, and 10^-10 otherwise.  It
%returns this in the field CTDout.eps.  It also returns the list of
%Overturns.  The field Overturns.keeper is set to one for each overturn
%which passsed all of the tests, and 0 for those that failed.  This routine
%stores NaN's in the keeper field, indicating that the second routine has
%not been called.  In this way, all overturns may be found with this
%routine, and EpsFromOverturns2 can be called multiple times with varied
%parameters, to see the effects of which overturns are kept and which are
%discarded.
%
%Inputs: 
%CTD has the following required fields:
%          lat: [1x75 double]
%          lon: [1x75 double]
%         year: [1x75 double]
%         yday: [1x75 double]
%            z: [2773x75 double]
%            T: [2773x75 double]
%            S: [2773x75 double]
%            H: [1x75 double]
%           id: [1x75 double]
%
%Notes: the fields T, S and C and D (optional) can have any name, specifiable in PP.  If D (density)
%is not specified, it is computed.  The fields lat, lon, year, yday, H and
%id are not required for this calculation but will not be passed to CTDout
%if they are not included, which will cause problems with plotting routines
%later.
%
%P should be in dbars.  S should be in psu.
%Data can be gridded or ungridded.
%
%PP A sample params struct is shown below.  Use the function
%DefaultOverturnPP to begin.
%
% PP.tvar='t';
% PP.svar='s';
% PP.cvar='c';
% PP.sgthvar='sgth'; %leave empty to compute it.
% PP.pvar='pgrid';
% PP.gridded='yes'; %yes if one pressure vector corresponds to all drops in structure
% PP.t_is_th='no'; %yes if T is potential temperature.
% PP.wh=1;
% PP.zmin=0;
% PP.zmax=400;
% PP.Startdown=10;
% PP.loud=0;
% PP.CUTOFF=1; %change later for run-length.
% PP.pref=0; %reference pressure in dbar
% PP.threshold=.01; %max Th for non-overturn.
% PP.plotit=1;
% 
%
%Outputs: CTDout is the de-Nan'ed, pressure-sorted version of CTD.  It has
%the additional fields D (if not specified earlier) and N2 (computed), and
%the sorted fields and Thorpe displacements as well.  The structure PP used
%to comptue it is stored as well.
%        Th: [2773x75 double]
%         D: [2773x75 double]
%         T: [2773x75 double]
%         S: [2773x75 double]
%       ThT: [2773x75 double]
%       ThC: [2773x75 double]
%        Ts: [2773x75 double]
%        Ds: [2773x75 double]
%        N2: [2773x75 double]
%         Z: [2773x75 double]
%         H: [1x75 double]
%       lat: [1x75 double]
%       lon: [1x75 double]
%        id: [1x75 double]
%      year: [1x75 double]
%        PP: [1x1 struct]
%
%Overturns: a list of structures, each of which has the fields shown in the
%sample below:

%            s: 1
%            f: 20
%           zs: 1
%           zf: 38.8000
%           N2: NaN
%         drho: 0.0019
%           dt: -0.0110
%           ds: 0.0050
%           Lt: 13.3860
%        Ltmax: -35.8000
%           Lp: 37.8000
%          eps: NaN
%          GKe: 0.5859
%     GKrunlen: 2.9439
%         Lt_t: 20.5030
%         Lt_c: 0
%           wh: 1
%       wh_out: 1
%          err: [1x1 struct]
%       keeper: 0
%
%The structure err is that output from GKerr, as discussed above - its
%fields are 
%          n: 2
%         dz: 1.9750
%       drho: 2.0000e-04
%          N: 0.0017
%         Lz: 3.9499
%       epsz: 5.2074e-08
%         Kz: 0.0035
%         Lr: 1.2705
%       epsr: 5.3873e-09
%         Kr: 3.5828e-04
%     ltsig4: 0.6660
%     ltsig3: 2.1061
%
%7/17/03, 
%MHA
%
%See also EpsFromOverturns2, ProfileOverturnPlot2,GKerr, runlength,
%OverturnCloseupPlot2




eval(['T=CTD.' PP.tvar ';']);
[m,n]=size(T);
if PP.wh > n | PP.wh==0 
    PP.wh=1:n;
end

if ~PP.loud
%    warning off MATLAB:polyfit:PolyNotUnique
    warning off % MATLAB:divideByZero
end

%Make an array for the output
CTDout.eps=zeros(m,length(PP.wh))*1e-10;

CTDout.Th=NaN*CTDout.eps;    
CTDout.D=NaN*CTDout.eps;
CTDout.T=NaN*CTDout.eps;
CTDout.D=NaN*CTDout.eps;
CTDout.S=NaN*CTDout.eps;
CTDout.ThT=CTDout.Th;
CTDout.ThC=CTDout.Th;
CTDout.Ts=CTDout.Th;
CTDout.Ds=CTDout.Th;
CTDout.N2=CTDout.Th;

if strcmp(PP.gridded,'yes')
    CTDout.Z=NaN*zeros(m,1);
else
    CTDout.Z=NaN*CTDout.eps;
end

Found=0; %Total overturns found

clear Overturn AllOverturns

for ii=1:length(PP.wh)
    wh=PP.wh(ii);
    %Get the fields out of the structure
    if strcmp(PP.gridded,'yes')        
        eval(['P=CTD.' PP.pvar ';']);
    else  
        eval(['P=CTD.' PP.pvar '(:,wh);']);
    end    
    %Get the desired range
    iz=find(P>PP.zmin & P<PP.zmax);
    P=P(iz);
    
    
    %truncate T and S to be the same
    eval(['T=CTD.' PP.tvar '(iz,wh);']);
    eval(['S=CTD.' PP.svar '(iz,wh);']);
    
    %Get cond too if requested
    if ~isempty(PP.cvar)
        eval(['C=CTD.' PP.cvar '(iz,wh);']);
    else 
        C=NaN*T;
    end
    
    %If T is potential temperature, copy it over.  If not, compute it. 
    if strcmp(PP.t_is_th,'yes')==1  
        TH=T;
    else
        %TH=sw_ptmp(S/1000,T,P/100,PP.pref/100); %recall that Mike rewrote all the toolbox routines 
        %to use MPa and s in c.u..
        TH=sw_ptmp(S,T,P,PP.pref); % not any more - 8/15
    end
    %If density has not been computed, compute it
    if ~isempty(PP.sgthvar)        
        eval(['D=CTD.' PP.sgthvar '(iz,wh);']);
    else
        % D=sw_dens(S/1000,TH,PP.pref/100);
        D=sw_dens(S,TH,PP.pref); % PSU and db
        if max(D)>1000
            D=D-1000; %sometimes the toolbox does not subtract 1000.  Do it here.
        end
    end
    
    %First remove all NaN's in the pressure
    inan=find(~isnan(P));
    if length(inan)<length(P)
        if PP.loud==1
            disp(['warning: removed ' num2str(length(P)-length(inan)) ' NaNs in pressure'])
        end
        D=D(inan);
        T=T(inan);
        S=S(inan);
        C=C(inan);
        P=P(inan);
    end
    
    
    
    %Now detect any pressure inversions and compute depth from pressure.
    %Also make the profile increasing in pressure (e.g. for SWIMS drops.)
    %Z=P;
    
    [Z,isp]=sort(P);
    T=T(isp);
    S=S(isp);
    D=D(isp);
    C=C(isp);
    
    %Now we have Z, T,S, D and C (C may be all NaN's.)
    
    %First remove all NaN's in the density from all variables
    inan=find(~isnan(D));
    if length(inan)<length(D)
        disp(['warning: removed ' num2str(length(D)-length(inan)) ' NaNs'])
        D=D(inan);
        T=T(inan);
        S=S(inan);
        Z=Z(inan);
        C=C(inan);
    end
    
    %At this point we could filter if we need to.
    %D=medfilt1(D,4);
    %Now Thorpe sort.  
    [Ds,is]=sort(D);
    [Ts,ist]=sort(T);
    %if temperature decreases with depth, reverse.
    if mean(diff(T))<0
        Ts=Ts(end:-1:1);
        ist=ist(end:-1:1);
    end
    ThT=Z(ist) - Z;
    
    %% NEW, June 2005:  add/remove density noise (DPW)
    np_typ = 'randn'; np_val = 0;
    if isfield(PP,'d_noise_proc')
        if iscell(PP.d_noise_proc)
            np_val = PP.d_noise_proc{1};
            np_typ = PP.d_noise_proc{2};
        else
            np_val = PP.d_noise_proc(1);
        end
    end
    % np_typ = 'randn' add random noise; 'sub' decrease flucts;
    %    'add' increase flucts;  np_val = sigma_theta noise fluct.
    if np_val>0
        switch np_typ
            case 'randn'
                D = D + (randn(size(D)) * np_val);
            case 'add'
                dD = D-Ds; % maximize flucts by specified value
                ix = find(dD<0); dD(ix) = dD(ix) - np_val;
                ix = find(dD>0); dD(ix) = dD(ix) + np_val;
                D = Ds + dD;
            case 'sub'
                dD = D-Ds; % minimize flucts, zero those smaller than spec
                ix = find(dD<0); dD(ix) = dD(ix) + np_val;
                iy = find(dD(ix)>0); 
                if ~isempty(iy)
                    dD(ix(iy)) = 0;
                end
                ix = find(dD>0); dD(ix) = dD(ix) - np_val;
                iy = find(dD(ix)<0); 
                if ~isempty(iy)
                    dD(ix(iy)) = 0;
                end
                D = Ds + dD;
        end
        % now, re-sort altered density profile:
        [Ds,is]=sort(D);
    end
    %% end of NEW stuff 
    % compute Thorpe displacements based on density:
    Th=Z(is) - Z;
    
    %repeat for cond if we have it
    [Cs,isc]=sort(C);
    %if cond decreases with depth, reverse.
    if mean(diff(C))<0
        Cs=Cs(end:-1:1);
        isc=isc(end:-1:1);
    end
    
    ThC=Z(isc) - Z;
    
    %Now get a list of overturns in this profile.
    Overturn=OverturnListFCN(Z,T,S,Th,D,Ds,PP);
    
    %Get depth interval for error calc
    dz=mean(diff(Z));
    
    %Fill in epsilon in the overturn locations.  
    %Here we enter the RMS Thorpe displacement from temperature
    %and conductivity too.
    for who=1:length(Overturn)
        whs=Overturn(who).s;
        whf=Overturn(who).f;
        Overturn(who).Lt_t=std(ThT(whs:whf));
        Overturn(who).Lt_c=std(ThC(whs:whf));
        Overturn(who).wh=wh;
        Overturn(who).wh_out=wh-PP.wh(1)+1;
        %error magnitudes
        dzOv = mean(abs(diff(Z(whs:whf))));
        err=GKerr(sqrt(Overturn(who).N2),0,PP.drho,dzOv,PP.n);
        Overturn(who).err=err;
        Overturn(who).keeper=NaN; %store a Nan since we have not yet decided whether it is a keeper
        %And fill in the dissipation - this should be done with 
        %the function EpsFromOverturns2.
        %CTDout.eps(whs:whf,wh-PP.wh(1)+1)=Overturn(who).eps;
    end
    
    %AllOverturns(Found+1:Found+length(Overturn))=Overturn;
    for c=1:length(Overturn)
        AllOverturns(Found+c)=Overturn(c);
    end
    Found=Found+length(Overturn);
    
    %Now fill in output structure
    
    %if we are not gridded, fill in all pressure vectors.  If not, only the
    %first
    if strcmp(PP.gridded,'yes')
        whp=1;
    else
        whp=wh-PP.wh(1)+1;
    end    
    
    CTDout.Z(1:length(Z),whp)=Z;
    %Fill in all the other variables
    CTDout.D(1:length(D),wh-PP.wh(1)+1)=D;
    CTDout.Ds(1:length(D),wh-PP.wh(1)+1)=Ds;
    CTDout.Ts(1:length(D),wh-PP.wh(1)+1)=Ts;
    CTDout.T(1:length(D),wh-PP.wh(1)+1)=T;
    CTDout.S(1:length(D),wh-PP.wh(1)+1)=S;
    CTDout.Th(1:length(D),wh-PP.wh(1)+1)=Th;
    CTDout.ThT(1:length(D),wh-PP.wh(1)+1)=ThT;
    CTDout.ThC(1:length(D),wh-PP.wh(1)+1)=ThC;
    
    %put in other fields into the output structure
    if isfield(CTD,'H')
        CTDout.H=CTD.H(PP.wh);
    end
    if isfield(CTD,'lat')
        CTDout.lat=CTD.lat(PP.wh);
        CTDout.lon=CTD.lon(PP.wh);
    end
    if isfield(CTD,'id')
        CTDout.id=CTD.id(PP.wh);
    end
    if isfield(CTD,'year')
        CTDout.year=CTD.year(PP.wh);
    end
    if isfield(CTD,'yday')
        CTDout.yday=CTD.yday(PP.wh);
    end

    CTDout.PP=PP;
    
    %dp=mean(diff(Z))*2;
    dp=10;
    if max(Z) < 10*dp
        dp=max(Z)/10;
    end
    % SKIP N^2 for now - DPW
%     [n2,p_n2]=nsqfcn(S/1000,T,Z/100,dp/100,dp/100);
%     CTDout.N2(1:length(Z),wh-PP.wh(1)+1)=interp1(p_n2*100,n2,Z);
    %A plot to examine the overturns
    if PP.plotit==1
        for who=1:length(Overturn);
            
            dzb=1;
            
            whs=Overturn(who).s-dzb;
            whf=Overturn(who).f+dzb;
            if whs<1
                whs=1;
            end
            if whf>length(D)
                whf=length(D);
            end
            
            figure(2)
            clf
            ax=MySubplot(.1,.3,0,.1,.1,0,4,1);
            axes(ax(1))
            plot(D(whs:whf),Z(whs:whf),Ds(whs:whf),Z(whs:whf))
            axis ij
            xlabel('D')
            title(['\epsilon= ' num2str(Overturn(who).eps) ', Lt= ' num2str(Overturn(who).Lt) ', GKe= ' num2str(Overturn(who).GKe)])
            
            axes(ax(2))
            plot(T(whs:whf),Z(whs:whf),Ts(whs:whf),Z(whs:whf))
            axis ij
            ytloff
            xlabel('T')
            axes(ax(3))
            plot(S(whs:whf),Z(whs:whf))
            axis ij
            ytloff
            xlabel('S')
            
            axes(ax(4))
            plot(Th(whs:whf),Z(whs:whf),ThC(whs:whf),Z(whs:whf),ThT(whs:whf),Z(whs:whf))
            axis ij
            xlabel('Thorpe')
            axes('position',[.72 .4 .24 .24])
            plot(S(whs:whf),T(whs:whf),S(whs+dzb:whf-dzb),T(whs+dzb:whf-dzb))
            axis square
            
            title(['Profile ' num2str(wh) ', overturn ' num2str(who) ' of ' num2str(length(Overturn))])
            
            pause
        end
    end
end

if ~exist('AllOverturns')
    AllOverturns=[];
end
