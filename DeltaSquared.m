function d2=DeltaSquared(f,r,T,rho)

W=2*pi/(3600*T);

a=r/((1-f).^(1/3));
c=r/((1-f).^(1/3))-f*r/((1-f).^(1/3));

d2=(EllPotTot(a,c,0,0,c,rho,W)-EllPotTot(a,c,a,0,0,rho,W)).^2;