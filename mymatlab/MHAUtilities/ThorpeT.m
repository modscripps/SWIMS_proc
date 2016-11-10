%function thorpeT=Thorpe(T,Tsort,DisplayP,DisplayZ);
%Use this function for temperature profiles, which decrease with depth.
%
%Return the Thorpe displacements (in bins, not meters) of a matrix T,
%whos columns represent temperature profiles, and its sorted matrix Tsort,
%sorted so that its high values are first, like temperature.
%(use 'Ts=sort (T); Ts=Ts(N:-1:1,:);').
%DisplayP and DisplayZ are optional parameters that allow for
%display of progress of the calculation.  Default is no notification.

function thorpeT=Thorpe(D,Dsort,DisplayP,DisplayZ);

%begin code here.  Densities are in matrix D, each col being one drop.
%Dsort=sort(D);
if (nargin < 4 )
	DisplayZ=0;
end
if (nargin < 3)
	DisplayP=0;
end

[m,n]=size(D);
thorpeT=zeros(m,n);
for c=1:n
	if (DisplayP ~= 0)
		DisplayProgress(c,DisplayP)
	end
	
	for d=1:m
		%if true temp < sorted temp look downwards
		if (DisplayZ ~= 0)
			DisplayProgress(d,DisplayZ)
		end
		if D(d,c)-Dsort(d,c) < -1e-10
			f=1;
			while D(d,c)-Dsort(d+f,c) < -1e-10
				f=f+1;
			end
			thorpeT(d,c)=-f;
		%if true temp > sorted temp look upwards
		elseif D(d,c)-Dsort(d,c) > 1e-10
			f=-1;
			while D(d,c)-Dsort(d+f,c) > 1e-10
				f=f-1;
			end
			thorpeT(d,c)=-f;
		else
			
		end	%if
	end
end

