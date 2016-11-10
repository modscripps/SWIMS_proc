function Overturn=OverturnListFCN(Z,T,S,displ,D,Ds,PP)
%Given profiles of sorted depth, T, S, Thorpe displacement, sgth and sorted
%sgth, find all overturns and return them in a list of structures.
%These should be sorted to remove nonmotonic pressure, and should have all
%NaN's removed.  See FindOverturnsFCN2 for details.
%
%MHA 6/30/03
%

%The following code is adapted from Jody's thorpe3_swims.m
dispsum=cumsum(displ);

Overturn=[];
G=9.80655;
rhoo=1026;

%
% Select samples with sums > threshold
i_overturn=find(dispsum>PP.threshold);
% ??? | (displ~=0 & dispsum<threshold));
%
% Terminate processing if there are no overturns
if isempty(i_overturn)
    Overturn=[];
    return;
end

%
% Select first points of an overturn as the last points
% after a gap in the indices of i_overturn.  This will not pick
% up the start of the uppermost overturn.
i_first=i_overturn(find(diff(i_overturn)>1)+1); 
%


% Select the last points in an overturn as the last points
% before a gap in the indices of i_overturn.  This will not pick
% up the bottom of an overturn that ends before the end of the
% record.
i_last=i_overturn(find(diff(i_overturn)>1));
%
% Find the start of the first overturn by checking that
% the first value of i_first should be less than the first
% value of i_last.
if isempty(i_first) | isempty(i_last)
    Overturn=[];
    return;
end;

if i_first(1)>1 & i_last(1) < i_first(1)
  i_first;
  i_overturn(1);
  i_first=[i_overturn(1); i_first];
end
%
% Find the end of the last overturn by checking if the last
% value of i_last is less than the last value of i_first.

if i_last(length(i_last)) < i_first(length(i_first))
  i_last=[i_last; i_overturn(length(i_overturn))];
end

% reject overturns that are below the threshold % DPW 5/2005
ix = [];
for i=1:length(i_last)
    if max( abs(dispsum( i_first(i):i_last(i) )) ) <= PP.cumsum_min
        ix = [ix i];
    end
end
i_last(ix) = []; i_first(ix) = [];

%Experimental - MHA. Find those pairs that are the same, and scoot them
%apart by one.
%isame=find(i_first==i_last & i_first ~= 1 & i_last ~= length(displ));
%i_first(isame)=i_first(isame)-1;
%i_last(isame)=i_last(isame)+1;

dz=mean(diff(Z));

PP.extrabins=max(fix(PP.dz_extra./dz),1);
%Also try moving things apart by a fixed distance.
isame=find(i_first > PP.extrabins & i_last < length(displ) - PP.extrabins);
i_first(isame)=i_first(isame)-PP.extrabins;
i_last(isame)=i_last(isame)+PP.extrabins;

% Determine the pressure bounds of the overturns and the number of
% overturns
pr_lb=Z(i_first); pr_ub=Z(i_last);
n_overturns=length(i_first);

% Get temperatures and salinities at the top and bottom of the
% overturns; lb & ub correspond to pressure limits, not temp
% or sal limits.
temp_lb=T(i_first); temp_ub=T(i_last);
sal_lb=S(i_first); sal_ub=S(i_last);

% Compute nsq across the overturns using resorted temp & sal
pr_avg=(pr_lb+pr_ub)/2;
%
% Obtain sorted temp & sal at pr_lb & pr_ub
tl=T(i_first); tu=T(i_last);
sl=S(i_first); su=S(i_last);
%
% Compute potential temp & potential density at pr_avg
%   for samples at tops of overturns (pr_lb)

% thl=sw_ptmp(sl/1000,tl,pr_lb/100,pr_avg/100); 
% sgthl=sw_pden(sl/1000,thl,pr_lb/100,pr_avg/100);
% %   for samples at bottoms of overturns (pr_ub)
% thu=sw_ptmp(su/1000,tu,pr_ub/100,pr_avg/100);
% sgthu=sw_pden(su/1000,thu,pr_ub/100,pr_avg/100);
% 
% now for PSU, db in seawater routines:
thl=sw_ptmp(sl,tl,pr_lb,pr_avg); 
sgthl=sw_pden(sl,thl,pr_lb,pr_avg);
%   for samples at bottoms of overturns (pr_ub)
thu=sw_ptmp(su,tu,pr_ub,pr_avg);
sgthu=sw_pden(su,thu,pr_ub,pr_avg);
%
% Use difference in sigma_thetas at pr_avg to compute nsq
%nsq=G^2*(sgthu-sgthl)./(1e6*(pr_ub-pr_lb));
%n2=g/rho*drho/dz - Jody's expression is right but for the wrong reasons
%and definitely wrong if P is in dbar.

%Signs... pr_ub is the higher pressure.  
nsq=G/rhoo*(sgthu-sgthl)./(pr_ub-pr_lb);

% Calculate the rms, mean, and largest displacements in each overturn.
rms_displ=NaN*ones(n_overturns,1); 
max_displ=NaN*ones(n_overturns,1);
meansigmat=NaN*ones(n_overturns,1);
for ii=1:n_overturns
  jj=find(Z>=pr_lb(ii) & Z<=pr_ub(ii));
  if ~isempty(jj)
    rms_displ(ii)=std(displ(jj));
    abs_displ=abs(displ(jj));
    max_displacement=max(abs_displ);
    imax=find(abs_displ==max_displacement);
    max_displ(ii)=displ(jj(imax(1)));
  else
    rms_displ(ii)=NaN; abs_displ=NaN;  max_displ(ii)=NaN; meansigmat(ii)=NaN;
  end;
end

%Put NaN's in where N^2 is negative
nsq(find(nsq<=0))=NaN;

for c=1:n_overturns
    Overturn(c).s=i_first(c);
    Overturn(c).f=i_last(c);
    Overturn(c).zs=pr_lb(c);
    Overturn(c).zf=pr_ub(c);
    Overturn(c).N2=nsq(c);
    Overturn(c).drho=sgthl(c)-sgthu(c);
    Overturn(c).dt=temp_lb(c)-temp_ub(c);
    Overturn(c).ds=sal_lb(c)-sal_ub(c);
    
    Overturn(c).Lt=rms_displ(c);
    Overturn(c).Ltmax=max_displ(c);
    Overturn(c).Lp=pr_ub(c)-pr_lb(c);
    Overturn(c).eps=0.64*Overturn(c).Lt.^2.*Overturn(c).N2.^1.5;
    Overturn(c).GKe=GKTestOverturn(Overturn(c),T,S,D,Ds);
    Overturn(c).GKrunlen=GKRunLengthTestOverturn(Overturn(c),Ds-D); %Compute run length of Thorpe fluxtuation not displacement
    
     
end

