%------------------------------------------------------------------------------
%         bootstd.m            09-16-92                Ren-Chieh Lien
%
%         Calculate mean, and standard errors of mean of time series
%
%
%         function [lbx,mx,ubx] = bootstd(x,n,prob,npts,lx,ux);
%
%         default npts = length(x);
%                 lx = -1.0e+20
%                 ux = 1.0e+20
%
%         ex.      [lbx,mx,ubx] = bootstd(x,512,95);
%
%         n times shaffling will be performed.
%
%         If # of x points is less than 10, NaN are assigned to lbx and ubx
%
%------------------------------------------------------------------------------
          function [lbx,mx,ubx] = bootstd(x,n,prob,npts,lx,ux);

          if (nargin <=5);
              ux = 1.0e+20;
          end
          if (nargin <= 4);
             lx = -1.0e+20;
          end
          bad = find(x>ux | x < lx | isnan(x));
          if (~isempty(bad)); x(bad) = []; end
          if (nargin <=3);
             npts = length(x);
          end
          
          npts = min([npts length(x)]);

          if (npts <= 5);
             lbx = NaN; mx = NaN; ubx=NaN;
          else
 
          minnpts = npts-0.01;
          k = floor(rand(npts*n,1)*minnpts+1);
          nx = reshape(x(k),npts,n);
   
          lprob=(100-prob)/2;
          uprob=lprob+prob;
          lprob = lprob*n/100;
          uprob = uprob*n/100;

          meannx = mean(nx);
          mnr = sort(meannx);
          estx=table1([(0:n+1)' [mnr(1);col(mnr);mnr(n)]], ...
             [lprob 0.5*n uprob]');
             lbx = estx(1);
             mx = estx(2);
             ubx = estx(3);
             mx = mean(mnr);
       end