function y = nansum(x);
% NANSUM Takes the sum of data excluding NaN's.
%
%       Y = nansum(X) returns a column vector Y, each element of which is
%       the average of the columns in X.  If X has NaN's they are excluded 
%       from the average. If the whole column is NaN, a NaN is still returned.

if isempty(x)
  y=NaN;
  return;
end;

bad = find(isnan(x));
x(bad)=0*bad;
y=sum(x);

