function [fh,fval]=HydrostaticStateExact3Ax2l(r1,r2,T,rho1,rho2,f_guess)


% [fh,fval]=fminsearch(@(f) DeltaSquared3Ax2l(f,r1,r2,T,rho1,rho2),f_guess,...
%     optimset('TolFun',1e-10));

A=[1 -1 0 0;
   0 0 1 -1;
   0 -1 0 1;
   -1 0 1 0];

b=[0 0 0 0];

DeltaSquared3Ax2l(f_guess,r1,r2,T,rho1,rho2)


 [fh,fval]=fminsearchcon(@(f) DeltaSquared3Ax2l(f,r1,r2,T,rho1,rho2),...
     f_guess,[0 0 0 0],[1 1 1 1],A,b,[],optimset('TolFun',1e-10,'MaxFunEvals',10000));