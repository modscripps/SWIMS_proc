function y = nanmean(x);
% NANMEAN Takes the mean of data excluding NaN's.
%
%       Y = nanmean(X) returns a column vector Y, each element of which is
%       the average of the columns in X.  If X has NaN's they are excluded 
%       from the average. If the whole column is NaN, a NaN is still returned.

% if a column vector simply return the vector...
[m,n] = size(x);
if m<=1
  y=x;
  return;
end;
good = ~isnan(x); % 0 for Nan, 1 otherwise.
if isempty(find(good))
  y=NaN;
  return;
end;
n_av = sum(good);  % find out how many good in each row.
bad_ind = find(~good);
if isempty(bad_ind)
  y=mean(x);return;
end;
x(bad_ind) = zeros(size(bad_ind));
y=NaN*x(1,:);
good_ave=find(n_av>0);
y(good_ave) = sum(x(:,good_ave))./n_av(good_ave);
bad_ind = find(n_av ==0);
y(bad_ind) = NaN*ones(size(bad_ind));
