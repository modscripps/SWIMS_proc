function GKerr=GKerr(N,plotit,drho,dz,n)
%function GKerr=GKerr(N,plotit,drho,dz,n)
%Return the G-K 1996 errors in displacement, dissipation rate and diapycnal
%diffusivity.  These are 
%
%Lz: the minimum detectable overturn given depth resolution dz.
%epsz: min detectable dissipation rate, from eps=0.64Lt^2N^3.
%Kz: same for K=0.2eps/N^2.
%Lr: the spurious overturn size that results from density errors drho.
%epsr, Kr: corresponding dissipation and diffusivity.
%ltsig3, 4: expected Thorpe scales for signals of K=1e-3 and 1e-4.
%Inputs: 
%N, buoyancy freq in rad/s
%plotit: 0 for no plot, 1 for Lt plot, 2 for eps, K plot
%drho: density precision, typical 1e-3 kg/m^3.
%dz: sample interval.
%n: # of points in vertical req'd to resolve an overturn.  GK use 5 (default).
%
%7/2, MHA
%
g=9.8;
rho=1026;

if nargin < 5
    n=5;
end
%Store input params
GKerr.n=n;
GKerr.dz=dz;
GKerr.drho=drho;
GKerr.N=N;

%Depth criteria
GKerr.Lz=n*dz*ones(size(N));
GKerr.epsz=0.64*n^2*dz^2.*N.^3;
GKerr.Kz=0.2*GKerr.epsz./N.^2;

%density criteria
GKerr.Lr=2*g/rho*drho./N.^2;
GKerr.epsr=0.64*GKerr.Lr.^2.*N.^3;
GKerr.Kr=0.2*GKerr.epsr./N.^2;

%signal strengths
K=1e-4;
GKerr.ltsig4=sqrt(1/.13*K./N);
GKerr.ltsig3=sqrt(1/.13*10*K./N);


if plotit==1
    
h=loglog(N.^2,GKerr.Lz,N.^2,GKerr.Lr,N.^2,GKerr.ltsig4,N.^2,GKerr.ltsig3);
legend(h,['depth res, \Delta z=' num2str(dz) ' m'],['noise, \Delta \rho = ' num2str(drho)],'sig,K=10^{-4}','sig,K=10^{-3}')
xlabel('N^2 / s^{-2}')
ylabel('L_T')
elseif plotit==2
subplot(211)    
h=loglog(N.^2,GKerr.epsz,N.^2,GKerr.epsr);
legend(h,['depth res, \Delta z=' num2str(dz) ' m'],['noise, \Delta \rho = ' num2str(drho)])
xlabel('N^2 / s^{-2}')
ylabel('\epsilon')
subplot(212)
h=loglog(N.^2,GKerr.Kz,N.^2,GKerr.Kr);
legend(h,['depth res, \Delta z=' num2str(dz) ' m'],['noise, \Delta \rho = ' num2str(drho)])
xlabel('N^2 / s^{-2}')
ylabel('K')

end