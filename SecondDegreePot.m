function U=SecondDegreePot(JS,mu,x,y,z,Rref,T)

J2=JS(1);
C22=JS(2);

W=2*pi/(3600*T);

[lambda,fi,r]=cart2sph(x,y,z);

U=mu./r;

rratio = (Rref./r).^2;

J2contrib = -mu./r.*(1.5*sin(fi).^2-0.5).*rratio.*J2./NormCoef(2,0);
C22contrib =-mu./r.*3.*((sin(fi).^2)-1).*cos(2.*lambda).*...
    rratio.*C22./NormCoef(2,2);

Urot=RotPot(x,y,W);

U=U+J2contrib+C22contrib+Urot;