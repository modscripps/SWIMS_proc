function [l,u]=CohConf(deg,perc,c2)
%Compute the upper and lower confidence limits on squared coherence.
%deg=degrees of freedom; perc = percent confidence level desired.
%From Percival's notes.
%Onyl works for 95% level.
q=1-(1-perc/100)^(2/(deg-2));

l=tanh(atanh(sqrt(c2))-1/(deg-2)-1.96/sqrt(deg-2));
l(find (l<0))=0;
u=tanh(atanh(sqrt(c2))-1/(deg-2)+1.96/sqrt(deg-2));
u(find (u>1))=1;
