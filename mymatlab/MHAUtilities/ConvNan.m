function varsm=ConvNan(var,tbins,logsmooth);

%THis code smooths the log over tbins*2+1 bins in time.  Each bin thus represents an average of the log
%of all profiles found within that time range.  If none are found, NaN is stored.  Note the number of
%profiles going into each varies, and therefore the std. deviation.
%#of bins on each side to smooth over.
%if logsmooth is nonzero, the mean will be computed of the log.  Ie 10^(mean(log)).

%tbins=12;
if ~exist('logsmooth')
   logsmooth=0;
end

[m,n]=size(var);
varsm=zeros(m,n)*NaN;
if logsmooth==0
   for c=1+tbins:n-tbins
      for d=1:m
         tmp=var(d,c-tbins:c+tbins);
         ind=find(~isnan(tmp));
         if isempty(ind)
            varsm(d,c)=NaN;
         else
            if length(ind)==1
               varsm(d,c)=tmp(ind);
            else
               varsm(d,c)=mean(tmp(ind));
            end
         end        
      end
   end
else
   for c=1+tbins:n-tbins
      for d=1:m
         tmp=var(d,c-tbins:c+tbins);
         ind=find(~isnan(tmp));
         if isempty(ind)
            varsm(d,c)=NaN;
         else
            if length(ind)==1
               varsm(d,c)=tmp(ind);
            else
               varsm(d,c)=10.^mean(log10(tmp(ind)));
            end
         end        
      end
   end
end
