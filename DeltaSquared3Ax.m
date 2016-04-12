function d2=DeltaSquared3Ax(f,r,T,rho)

fp=f(1);
fq=f(2);

W=2*pi/(3600*T);

a=r/(((fp-1)*(fq-1))^(1/3));
b=(r - fq*r)/(((fp-1)*(fq-1))^(1/3));
c=(r - fp*r)/(((fp-1)*(fq-1))^(1/3));

d2=(EllPotTot3Ax(a,b,c,0,b,0,rho,W)-EllPotTot3Ax(a,b,c,a,0,0,rho,W)).^2+...
   (EllPotTot3Ax(a,b,c,0,0,c,rho,W)-EllPotTot3Ax(a,b,c,a,0,0,rho,W)).^2;

