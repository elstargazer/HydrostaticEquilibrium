function U=EllPotTot3Ax2l(a1,b1,c1,a2,b2,c2,x,y,z,rho1,rho2,W)

% Outer layer potential
% [a1 b1 c1]
% 
% [a2 b2 c2]

U=Ell3Pot(a1,b1,c1,x,y,z,rho1);

% Inner layer potential

U=U+Ell3Pot(a2,b2,c2,x,y,z,rho2-rho1);

% Rotational potential

U=U+RotPot(x,y,W);

% Tidal potential

U=U+TidalPot(x,y,z,W);


