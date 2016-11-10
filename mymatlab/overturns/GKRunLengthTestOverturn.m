function rmslen=GKTest(Overturn,Th)
%function runlength=GKTest(Overturn,Th)
%Compute the RMS run length of Thorpe fluctuation over an overturn.
s=Overturn.s;
f=Overturn.f;

data=Th(s:f);

[rmslen,locs]=runlength(data);