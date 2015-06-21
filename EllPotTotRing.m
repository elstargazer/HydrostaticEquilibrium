function U=EllPotTotRing(a,c,x,y,z,rho,W,Rr,Mr)

U=RotPot(x,y,W)+EllPot(a,c,x,y,z,rho)+RingPot(x,y,z,Rr,Mr);