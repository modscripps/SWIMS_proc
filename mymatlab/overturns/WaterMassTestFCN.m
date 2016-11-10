function [Dc,Starts, Ends, GoodStarts, GoodEnds,Regions, GoodRegions, RegionsCorrected, GKe]=WaterMassTestFCN(Z,T,S,D,Ds,PP)
%function [Dc,Starts, Ends, Regions, GoodRegions, RegionsCorrected]=WaterMassTestFCN(Z,T,S,D,Ds,PP)
%An updated version of the water-mass test.

Ndown=length(D);

Startdown=PP.Startdown;

Regions=0;
GoodRegions=0;
RegionsCorrected=0;
Starts=[]; Ends=[];
GoodStarts=[];GoodEnds=[];
GKe=[];
%We now have Z,T,S,D and Ds.

Dc=D;

c=PP.Startdown;
while c < Ndown
    %Search for start of reordering region
    while D(c) == Ds(c) & c < Ndown
        c=c+1;
    end
    
    %Don't go through rigamarole if we've only found
    %end of profile
    if c ~= Ndown
                
        %c now indexes the first point where Ds != D, or Ndown
        s=c-1;
        
        %so s indexes the last equal point
        if PP.loud==1
            disp 'Found start of reordering region at s='
            s
        end
        
        %Search for the end of the region
        while Ds(c) < max(D(s:c)) & c < Ndown
            c=c+1;
        end
        
        %c indexes the last unequal point so add one!
        c=c+1;
        
        %c now indexes first equal point past entire region, or Ndown if end was reached first
        f=min(c,Ndown);
        
        if PP.loud==1
            disp 'Found end of reordering region at f='
            f
        end
        
        %Now we must have the region start before Ndown - 1 since
        %we subtracted one.
        if s < Ndown - 1
            %1/23/97 change.  If reordering region extends past Ndown then still examine it.	
            Regions=Regions+1;
            e=0;
            TestRegion2002;
            if PP.loud == 1
                e
            end
            
            if e==0 | e > .5  | f - s < PP.CUTOFF
                %region rejected - replace with the sorted profile in the range.
                Dc(s:f)=Ds(s:f);
                RegionsCorrected=RegionsCorrected+1;
                Starts(Regions)=s;
                Ends(Regions)=f;
                GKe(Regions)=e;
            else
                
                GoodRegions=GoodRegions+1;
                GoodStarts(GoodRegions)=s;
                GoodEnds(GoodRegions)=f;
                Starts(GoodRegions)=s;
                Ends(GoodRegions)=f;
            end
        end
    end
end		%while c < Ndown

if PP.loud==1
    disp 'Regions Found'
    Regions
    disp 'Regions Corrected'
    RegionsCorrected
end

%Compute Thorpe
%[Dcs,is]=sort(Dc);
%Thc=Z(is) - Z;
