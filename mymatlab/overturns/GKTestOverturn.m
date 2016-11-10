function e=GKTest(Overturn,T,S,D,Ds)
%This is a newer version of the GKtest routine.  I replaced norm with std,
%which makes no difference, but more dignificantly it is important only to
%use data in the reordering region, not the outsides as in the N
%computation.  
%7/3/03
%MHA


warning off %MATLAB:polyfit:RepeatedPointsOrRescale
s1=Overturn.s;
f1=Overturn.f;

%Since we allow some extra on the ends when computing N, only use the
%reording region here.
e=10; % will cause to be rejected - DPW
ind1=find(Ds(s1:f1)~=D(s1:f1));
if length(ind1)>1
    s=s1+ind1(1);
    f=s1+ind1(end);
    
    %Solve Gm=d where G=(ones S's) m=(a b)' and d=densities
    ps=polyfit(S(s:f),D(s:f),1);
    pt=polyfit(T(s:f),D(s:f),1);
    
    rhos=polyval(ps,S(s:f));
    rhot=polyval(pt,T(s:f));
    
    rmsth=std(D(s:f)-Ds(s:f));
    
    es=std(rhos-D(s:f))./rmsth;
    et=std(rhot-D(s:f))./rmsth;
    e=max(es,et);
end
