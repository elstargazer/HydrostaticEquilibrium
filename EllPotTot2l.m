function U=EllPotTot2l(a1,c1,a2,c2,x,y,z,rho1,rho2,W)


% Outer layer potential



U=EllPot(a1,c1,x,y,z,rho1);


% Inner layer potential

U=U+EllPot(a2,c2,x,y,z,rho2-rho1);

% Rotational potential

U=U+RotPot(x,y,W);