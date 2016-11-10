function delp=PhaseConf(deg,perc,c2)
%Return the error bars on phase, in degrees, 
%from Dave's code. From Jenkins+Watts 1968 p 434 or Otnes+Enochson 1978 p 381
x1=fv95(2,deg-2,1-perc/100);

x2=2/(deg-2);

x3=c2;
x4=zeros(size(x3));
ind=find(x3>1e-20);
x4(ind)=sqrt(x2.*x1.*(1-x3(ind))./x3(ind));

delp=zeros(size(x4))+180;

ind=find(x4<1);

delp(ind)=abs(asin(x4(ind)))/pi*180;

function f=fv95(n1,n2,a1)

b=n1/2;
a=n2/2;
y=nv95(1-a1);
lam=(y.^2-3)/6;
h=2/(1/(2.*a-1)+1/(2.*b-1));

w=y.*sqrt(h+lam)./h-(1./(2.*b-1)-1./(2.*a-1)).*(lam+0.8333333333333333-2/3./h);

f=exp(2.*w);

function c=nv95(p)

c=p;
ind=find(1-p<p);
c(ind)=1-p(ind);
c=sqrt(-2.*log(c));

c=c-(2.30753+0.27061.*c)./(1+.99229.*c+0.04481.*c.^2);
c=c.*(p-.5)./abs(p-.5);