%function d=VarTrack2(var,depth,isopycdepth)
%
%VarTrack2.m
%Matthew Alford
%Modified from VarTrack2.m retrieved from 12/9/96 backup,
%on 1/27/97 when the original mysteriously disappeared.
%
%Track a variable into isopycnal coordinates by inputting 
%the depths of each isopycnal, previously tracked with VarTrack,
%the depth vector reflecting the depth of the Eulerian measurements,
%and the eulerian data.  Returns the semilagrangian variable. 
%
%INPUT var = [nbins x nrecords],depth=[nbins x 1], 
%isopycdepth = [Num_Isopycnals x nrecords]
%
%OUTPUT d = [Num_Isopycnals x nrecords]
%
%example: Ri_Iso=VarTrack2(Ri,z,depth);
%
function d=UnTrack(var,depth,isopycdepth)
[m,n]=size(var);
%m is the number of bins of var, n is the number of records.
N=length(depth);

%Modified for compatibility with NaN's
ind=find(isnan(var));
var(ind)=-9999;

%Also do some error checking on isopycdepth: find any < 0 values.
ind=find(isopycdepth<0);
isopycdepth(ind)=-9999;
ind=find(isnan(isopycdepth));
isopycdepth(ind)=-9999;


d=zeros(N,n);
for b=1:n
  d(:,b)=interp1(isopycdepth(:,b),var(:,b),depth,'linear');
end
  		
%One more change: now fill in NaN's wherever there were NaN's in isopycdepth
%d(find(d<0)) = NaN;