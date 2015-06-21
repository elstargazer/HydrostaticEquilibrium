function U=EllPotTot3Ax(a,b,c,x,y,z,rho,W)

 U=RotPot(x,y,W)+Ell3Pot(a,b,c,x,y,z,rho)+TidalPot(x,y,z,W);

% U=RotPot(x,y,W)+Ell3Pot4approx(a,b,c,x,y,z,rho);