function PP=DefaultOverturnPP;
%function PP=DefaultOverturnPP;
%Return a list of default parameters for a Thorpe Scale calculation.
%Names of fields
PP.tvar='T'; 
PP.svar='S';
PP.cvar='';
PP.sgthvar=''; %leave empty to compute it.
PP.pvar='z';
PP.gridded='no'; %yes if one pressure vector corresponds to all drops in structure
PP.t_is_th='no'; %yes if T is potential temperature.
PP.wh=1; %The range of drops to do.
PP.zmin=0; %A depth range to resrict calculation to.
PP.zmax=5000;

PP.loud=0; %certain diagnostics will be output to screen if this is 1.
PP.CUTOFF=7; %Galbraith-Kelley Run-len criterion.
PP.pref=0.01; %reference pressure in dbar

PP.threshold=1; % min cumsum(Th) to delineate a potential overturn.
PP.cumsum_min=PP.threshold; % max cumsum(Th) for non overturn (reject potential ones)

PP.plotit=0; %Some plots will be output if this is 1.

PP.THRESH_t=0; %The fraction of the density Thorpe scale that the temperature Thorpe scale must be to count as an overturn.
PP.THRESH_c=0; %same as above, for conductivity.
PP.THRESH_GK=1.2; %Put to a huge # to disable GK checks.  Formally, GK say this should be 0.5.

%Error parameters - density precision and the # of samples needed to
%resolve an overturn
PP.drho=1e-3; %estaimte of density precision
PP.n=2; %estimate of # of samples req'd to resolve an overturn - GK say 5.

PP.dz_extra=.1; %distance in m to include at top and bottom of each overturn in computing N

PP.THRESH_r=1; %The factors by which the observed overturn scale must exceed depth and density overturn scales.
PP.THRESH_z=1; %
