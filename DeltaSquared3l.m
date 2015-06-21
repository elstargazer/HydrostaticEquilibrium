function d2=DeltaSquared3l(f,r1,r2,r3,T,rho1,rho2,rho3)

W=2*pi/(3600*T);

a1=r1/((1-f(1)).^(1/3));
c1=r1/((1-f(1)).^(1/3))-f(1)*r1/((1-f(1)).^(1/3));

a2=r2/((1-f(2)).^(1/3));
c2=r2/((1-f(2)).^(1/3))-f(2)*r2/((1-f(2)).^(1/3));

a3=r3/((1-f(3)).^(1/3));
c3=r3/((1-f(3)).^(1/3))-f(3)*r3/((1-f(3)).^(1/3));

% 1 outer
% 2 inner

d2=(EllPotTot3l(a1,c1,a2,c2,a3,c3,0,0,c1,rho1,rho2,rho3,W)-EllPotTot3l(a1,c1,a2,c2,a3,c3,a1,0,0,rho1,rho2,rho3,W)).^2+...
   (EllPotTot3l(a1,c1,a2,c2,a3,c3,0,0,c2,rho1,rho2,rho3,W)-EllPotTot3l(a1,c1,a2,c2,a3,c3,a2,0,0,rho1,rho2,rho3,W)).^2+...
   (EllPotTot3l(a1,c1,a2,c2,a3,c3,0,0,c3,rho1,rho2,rho3,W)-EllPotTot3l(a1,c1,a2,c2,a3,c3,a3,0,0,rho1,rho2,rho3,W)).^2;