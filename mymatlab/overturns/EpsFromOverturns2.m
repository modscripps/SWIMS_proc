function [Overturn,CTDout]=EpsFromOverturns2(Overturn, CTDout)
%See also FindoverturnsFCN2.  This function takes the list of overturns and
%the CTDout structure from that function, and, using the parameters in
%CTDout.PP, computes dissipation by storing the dissipation for each
%overturn in the profile for the overturn's range spanned.  It stores a 1
%in Overturn.keeper for each one kept, and 0 for the ones thrown out.
%
%Possible tests are: temperature inversion too, Galbraith-Kelley watermass
%and run-length tests, RMS Thorpe displacement greater than the expected
%noise from pressure and density.
%
%7/17/03
%MHA
%

PP=CTDout.PP;

CTDout.eps=zeros(size(CTDout.T))*1e-10;
% DPW 8/03, added krho, L_Th for each valid overturn
CTDout.krho=zeros(size(CTDout.T))*5e-6;
CTDout.L_Th=zeros(size(CTDout.T))*0;
for who=1:length(Overturn)
    %if we're gridded, of course just use the single pressure profile.
    %Otherwise, use the proper one.
    if strcmp(PP.gridded,'no')
        %iprof=Overturn(who).wh; %whoops - use wh_out here and below
        iprof=Overturn(who).wh_out;
    else
        iprof=1;
    end
    io=find(CTDout.Z(:,iprof)>Overturn(who).zs & CTDout.Z(:,iprof) < Overturn(who).zf);
    
    
    if Overturn(who).Lt_t >= PP.THRESH_t * Overturn(who).Lt & ...
            Overturn(who).Lt_c >= PP.THRESH_c * Overturn(who).Lt &  ...      %temperature inversion too
            Overturn(who).GKe <= PP.THRESH_GK &  ...                        %GK OK
            Overturn(who).Lt > PP.THRESH_r*Overturn(who).err.Lr & ...       %Greater than the detectable in terms of rho
            Overturn(who).Lt > PP.THRESH_z*Overturn(who).err.Lz & ...       % "         "       "           "       z    
            Overturn(who).GKrunlen > PP.CUTOFF  & ...                       %G-K run length test  
            ~isnan(Overturn(who).N2)
        
        Overturn(who).keeper=1;
        CTDout.eps(io,Overturn(who).wh_out)=Overturn(who).eps;
        CTDout.krho(io,Overturn(who).wh_out)=0.2*Overturn(who).eps./Overturn(who).N2;
        CTDout.L_Th(io,Overturn(who).wh_out)=Overturn(who).Lt;
    else
        Overturn(who).keeper=0;
    end
    
end