function [ydaysm,prsm,varsm]=SmoothVar(varstr,BinsPerDay,nv,logsm)%function [ydaysm,prsm,varsm]=SmoothVar(varstr,BinsPerDay,nv,logsm)%%Load a gridded BS98 variable and smooth it.%Var can be Eps,Obs,Sal, Theta, Sigmat or Krho. d_Iso works too.%This computes the average in centered bins of the variable.%NaN's are not averaged in.  Output is nv by 14*BinsPerDay.%If logsm is one, the logarithm is smoothed instead of arithmetic mean.%default is 0.if ~exist('logsm')	logsm=0;end%Compute krho is that's what we wantif strcmp(varstr,'Krho')==1   %Now we get krho%   load 'E:\Data\mmp\bs98\gridded\Eps'%   load 'E:\Data\mmp\bs98\gridded\N2'   load 'Shortboard:Data:mmp:bs98:gridded:Eps'   load 'Shortboard:Data:mmp:bs98:gridded:N2'      ind=find(N2>=0);   var=Eps.*NaN;   var(ind)=0.2*Eps(ind)./N2(ind);else %otherwise we just want one variable    %   loadstr=['load E:\Data\mmp\bs98\gridded\' varstr];	loadstr=['load Shortboard:Data:mmp:bs98:gridded:' varstr];	eval(loadstr)      eval(['var= ' varstr ';'])      if strcmp(varstr,'N2')==1      ind=find(var<0);      var(ind)=NaN;   endend%Now grid data.firstday=294;lastday=307;days=lastday-firstday+1;%# of bins to average overvbins=fix(300/nv);varsm=zeros(nv,BinsPerDay*days)*NaN;prm=zeros(vbins,1);for c=1:BinsPerDay*days   DisplayProgress(c,20)   indt=find(yday > firstday+ (c-1)/BinsPerDay & yday <= firstday + c/BinsPerDay);   if isempty(indt)      %no profiles in this time slot      varsm(:,c)=NaN;   else      for d=1:nv			if logsm==0         	tmp=var((d-1)*vbins+1:d*vbins,indt);         else         	tmp=log10(var((d-1)*vbins+1:d*vbins,indt));			end			ind=find(~isnan(tmp));         if isempty(ind)            varsm(d,c)=NaN;         else            if length(ind)==1               varsm(d,c)=tmp(ind);            else               varsm(d,c)=mean(tmp(ind));            end         end               end   endendif logsm==1	varsm=10.^varsm;endfor d=1:nv   prsm(d)=mean(pr((d-1)*vbins+1:d*vbins));endydaysm=firstday+((1:BinsPerDay*days)-0.5)/BinsPerDay;