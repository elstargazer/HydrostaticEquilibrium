function U=EllPotTot3AxRing(a,b,c,x,y,z,rho,W)

Rr=900000;
Mr=4e29;
1
 U=RotPot(x,y,W)+Ell3Pot(a,b,c,x,y,z,rho)+RingPot(x,y,z,Rr,Mr);

% U=RotPot(x,y,W)+Ell3Pot4approx(a,b,c,x,y,z,rho);

