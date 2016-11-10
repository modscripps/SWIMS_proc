%function thorpe=Thorpe(D,Dsort);

function thorpe=Thorpe(D,Dsort);

%Return the Thorpe displacements (in bins, not meters) of a matrix D,
%whos columns represent density profiles, and its sorted matrix Dsort
%(use 'sort').

%begin code here.  Densities are in matrix D, each col being one drop.
%Dsort=sort(D);
[m,n]=size(D);
thorpe=zeros(m,n);
for c=1:n
	for d=1:m
		%if true dens > sorted dens look downwards
		if D(d,c)-Dsort(d,c) > 1e-10
			f=1;
			while D(d,c)-Dsort(d+f,c) > 1e-10
				f=f+1;
			end
			thorpe(d,c)=-f;
		%if true dens < sorted dens look upwards
		elseif D(d,c)-Dsort(d,c) < -1e-10
			f=-1;
			while D(d,c)-Dsort(d+f,c) < -1e-10
				f=f-1;
			end
			thorpe(d,c)=-f;
		else
			
		end	%if
	end
end

