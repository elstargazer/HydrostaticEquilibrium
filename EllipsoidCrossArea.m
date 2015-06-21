function [area]=EllipsoidCrossArea(r,fq,fp,fi,ph)

th=pi/2-fi;

a=r/(((fp-1)*(fq-1))^(1/3));
b=(r - fq*r)/(((fp-1)*(fq-1))^(1/3));
c=(r - fp*r)/(((fp-1)*(fq-1))^(1/3));


a2=a*a;
b2=b*b;
c2=c*c;


cth=cos(th);
sth=sin(th);
sph=sin(ph);
cph=cos(ph);

sth2=sth.*sth;
cth2=cth.*cth;
sph2=sph.*sph;
cph2=cph.*cph;
aa=cth2.*(cph2./a2+sph2./b2)+sth2./c2;
twohh=2.*cth.*sth.*cph*(1./b2-1./a2);
bb=sph2./a2+cph2./b2;
ps=0.5*atan2(twohh,aa-bb);
sps=sin(ps);
cps=cos(ps);
aaa=cps.*(aa.*cps+twohh.*sps)+bb.*sps.*sps;
bbb=sps.*(aa.*sps-twohh.*cps)+bb.*cps.*cps;
semax1=1./sqrt(aaa);
semax2=1./sqrt(bbb);

area=pi.*semax1.*semax2; 

