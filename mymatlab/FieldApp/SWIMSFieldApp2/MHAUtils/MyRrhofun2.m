function [Rrho,angle]=MyRrhofun2(s,th,p,zsm);%Rrhofun2: [Rrho,angle]=MyRrhofun2(s,th,p) - density ratio as angle%Smooth over interval zsm first.%Note this differs from Mike's: my definition Tu(me) = Tu(mike) - 45 degs.% Rrho is returned as well as Tu.  Note Rrho is alphaTz/betaSz, not%the other way around.  This is tested and works now.%s in psu, and p in MPa.%sample interval in mzs=100*nanmean(diff(p));%Make a filter and smooth[b,a]=MHAbutter(zs,zsm);ind=~isnan(s);ss=zeros(size(s));ss(ind)=filtfilt(b,a,s(ind));ths=zeros(size(th));ths(ind)=filtfilt(b,a,th(ind));len=length(p);alpha=sw_alpha(ss,ths,100*p);beta=sw_beta(ss,ths,100*p);%z is positive upwards.  So when s is increasing with depth dsdz < 0.%t decreasing with depth --> dtdz > 0.  So dp=diffs(p);indg=find(dp~=0);dsdz=zeros(size(p))*NaN;dthdz=dsdz;dsdz(indg) = -diffs(ss(indg))./dp(indg);dthdz(indg) = -diffs(ths(indg))./dp(indg);%    bs = -beta.*dsdz;%    at = alpha.*dthdz;%	angle=atan((at - bs)./(at+bs));%12/30/99 changes below.  Turner angle was mapped improperly%when it was less than -pi/2.bs = beta.*dsdz;at = alpha.*dthdz;%Now make all vectors the size of the inputs.[m,n]=size(s);Rrho=zeros(m,n)*NaN;indg=find(bs~=0);%Rrho(1:m,:)Rrho(indg)=at(indg)./bs(indg);angle=atan(-Rrho)-pi/4;	ipi=find(angle<-pi/2);angle(ipi)=angle(ipi)+pi;%angle=atan2(at+bs,at-bs);