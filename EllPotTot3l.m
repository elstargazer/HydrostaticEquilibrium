function U=EllPotTot3l(a1,c1,a2,c2,a3,c3,x,y,z,rho1,rho2,rho3,W)


% 1 Outer layer potential



U=EllPot(a1,c1,x,y,z,rho1);


% middle layer potential

U=U+EllPot(a2,c2,x,y,z,rho2-rho1);

% inner layer potential

U=U+EllPot(a3,c3,x,y,z,rho3-rho2);

% Rotational potential

U=U+RotPot(x,y,W);