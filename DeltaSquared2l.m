function d2=DeltaSquared2l(f,r1,r2,T,rho1,rho2)

W=2*pi/(3600*T);

a1=r1/((1-f(1)).^(1/3));
c1=r1/((1-f(1)).^(1/3))-f(1)*r1/((1-f(1)).^(1/3));

a2=r2/((1-f(2)).^(1/3));
c2=r2/((1-f(2)).^(1/3))-f(2)*r2/((1-f(2)).^(1/3));

% 1 outer
% 2 inner

d2=(EllPotTot2l(a1,c1,a2,c2,0,0,c1,rho1,rho2,W)-EllPotTot2l(a1,c1,a2,c2,a1,0,0,rho1,rho2,W)).^2+...
  (EllPotTot2l(a1,c1,a2,c2,0,0,c2,rho1,rho2,W)-EllPotTot2l(a1,c1,a2,c2,a2,0,0,rho1,rho2,W)).^2;