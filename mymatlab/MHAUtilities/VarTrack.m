%function d=VarTrack(var,dens,isopyc)
%
%VarTrack.m
%Matthew Alford
%Modified from VarTrack.m retrieved from 12/9/96 backup,
%on 1/27/97 when the original mysteriously disappeared.
%
%INPUT var = [nbins x nrecords],dens = [nbins x nrecords],
%isopyc = [1 x Num_Isopycnals]
%
%OUTPUT d = [Num_Isopycnals x nrecords]
%
%Take a variable to track, and the density at each bin,
%and return a list of that variable at each of a list of input
%isopycnals.
%
%example d_Iso=VarTrack(Z,D,newpycs);
%
function d=VarTrack(var,dens,isopyc)

[m,n]=size(var);
%m is the number of bins of var, n is the number of records.
%dens had better be the same size.
[M,N]=size(isopyc);

%Modified for compatibility with NaN's
ind=find(isnan(var));
var(ind)=-9999;
ind=find(isnan(dens));
dens(ind)=-9999;

d=zeros(N,n);
for b=1:n
	ind=1;
	for c=1:N
		ind=ind-10;
		
		if (ind < 1)
			ind = 1;
		end

		while dens(ind,b) < isopyc(c) &  ind < m 
			ind=ind+1;
		end
		ind=ind-1;
	
		%ind should now index the density just smaller than the isopycnal we
		%seek.
		
		%if the sought isopyc is too big return 500, if too small return 0.
		if ind == m - 1 
			d(c,b)=5000;
		elseif ind == 0
			d(c,b)=0;
		else
			diff1=(dens(ind+1,b) - dens(ind,b));
			if diff1 == 0.0 
				diff1=1e-10;
			end
			fract=( isopyc(c)-dens(ind,b)) / diff1;
			d(c,b)=var(ind,b)+fract*(var(ind+1,b)-var(ind,b));
		end
	end
end
		
