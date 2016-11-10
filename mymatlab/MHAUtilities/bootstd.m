%------------------------------------------------------------------------------
%         bootstd.m            09-16-92                Ren-Chieh Lien
%
%         Calculate mean, and standard errors of a time series
%
%
%         function [lbx,mx,ubx] = bootstd(x,n,prob);
%
%         ex.      [lbx,mx,ubx] = bootstd(x,512,95);
%
%         n times shaffling will be performed.
%
%         If # of x points is less than 10, NaN are assigned to lbx and ubx
%
%------------------------------------------------------------------------------
          function [lbx,mx,ubx] = bootstd(x,n,prob);

          bad = find(isnan(x));
          x(bad) = [];
          npts = length(x);
          if ( n >= 500); n = 500; end
          if (npts <= 10);
             lbx = NaN; mx = NaN; ubx=NaN;
          else
             minnpts = npts-0.01;
             k = floor(rand(npts*n,1)*minnpts+1);
             nx = reshape(x(k),npts,n);
   
             lprob=(100-prob)/2;
             uprob=lprob+prob;
             lprob = lprob*npts/100;
             uprob = uprob*npts/100;

             nr = sort(nx);
             mnr = mean(nr')';
             estx=interp1((0:npts+1)',[mnr(1);mnr(:);mnr(npts)], ...
                   [lprob 0.5*npts uprob]');

             lbx = estx(1);
             mx = estx(2);
             ubx = estx(3);
       end