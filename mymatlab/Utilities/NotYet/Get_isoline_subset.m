function [xsub, ysub] = Get_isoline_subset(x, y, xlim, ylim)
% function [xsub, ysub] = Get_isoline_subset(x, y, xlim, ylim);
% Use to retrieve (x,y) values for an isoline, retaining NaNs
%   to allow proper plotting, EG, of isobaths, shorelines with
%   separate closed contours (islands, etc).
% x, y = entire set (with NaNs at gaps),
% xlim, ylim = range specifications,
% xsub, ysub = desired subset, with gaps retained
% DPW, apr 2002

iok = find( (x>=xlim(1)&x<=xlim(2) & y>=ylim(1)&y<=ylim(2)) | isnan(x+y));
xsub = x(iok); ysub = y(iok);

% EG: ib = find([IsoBath(:).depth] == dep);
% then dep-isobath is at (IsoBath(ib).lon,IsoBath(ib).lat),
% use this function to get subset for plot limits