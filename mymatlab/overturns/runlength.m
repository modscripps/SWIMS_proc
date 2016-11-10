function [rmslen,locs]=runlength(data)
%Compute the locations, locs, of zero crossings in a data series data.  
%7/03/03 MHA

%
%Get the sign of the first point
first=min(find(data~=0));

if isempty(first) | first > length(data)-1
    rmslen=NaN;locs=[];return
end

nowsgn=sign(data(first));
locs(1)=first;
crossings=1;
for c=first+1:length(data)
    if sign(data(c)) ~= nowsgn & sign(data(c)) ~=0
        crossings=crossings+1;
        locs(crossings)=c;
     
        nowsgn=sign(data(c));
    end
end

%Now if the last point was not a zero crossing, add on the last point to
%get the last 'run'.
% if locs(end)~=length(data)
%     crossings=crossings+1;
%     locs(crossings)=length(data);
% end

rmslen=sqrt(mean(diff(locs).^2));