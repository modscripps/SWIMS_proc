%function d=VarTrack3(var,depth,isopycdepth)
%
%VarTrack3.m
%Matthew Alford
%Modified from VarTrack2.m
%
%Track a variable into isopycnal coordinates by inputting 
%the depths of each isopycnal, previously tracked with VarTrack,
%the depth vector reflecting the depth of the Eulerian measurements,
%and the eulerian data.  Returns the semilagrangian variable. 
%Differs from VarTrack2 in that the depth changes with time;
%ie, it is a matrix not a vector, as in CTD data.
%
%INPUT var = [nbins x nrecords],depth=[nbins x nrecords], 
%isopycdepth = [Num_Isopycnals x nrecords]
%
%OUTPUT d = [Num_Isopycnals x nrecords]
%
%example: Cx_Iso=VarTrack3(Cx,Nb,depth);
%
function d=VarTrack2(var,depth,isopycdepth)

[m,n]=size(var);
%m is the number of bins of var, n is the number of records.
%depth had better be the same number of bins, but a vector.
%Also isopycdepth is a matrix not a vector-- and each column is one time.
%So we replace N by M below.
[M,N]=size(isopycdepth);
d=zeros(N,n);
for b=1:n
	ind=1;
	for c=1:M
		ind=min(find(depth(:,b) > isopycdepth(c,b)))-1;
		%ind should now index the depth just smaller than the isopycnal we
		%seek.
		
		%if the sought isopyc is too big return 500, if too small return 0.
		if isempty(ind)	%ind == [] ; changed for version 5 3/18/97
			d(c,b)=5000;
		elseif ind == 0
			d(c,b)=0;
		else
			diff1=(depth(ind+1,b) - depth(ind,b));
			if diff1 == 0.0 
				diff1=1e-10;
			end

			fract=( isopycdepth(c,b)-depth(ind,b)) / diff1;
			d(c,b)=var(ind,b)+fract*(var(ind+1,b)-var(ind,b));
		end
	end
end
		
