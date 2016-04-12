function [fh,fval]=HydrostaticStateExact2lTidal(r1,r2,T,rho1,rho2,f0)

A=[-1 1 0 0;
    0 0 -1 1;
    0 -1 0 1;
    -1 0 1 0];

b=[0 0 0 0];

[fh,fval]=fminsearchcon(@(f) DeltaSquared2lTidal(f,r1,r2,T,rho1,rho2),...
    f0,[0 0 0 0],[1 1 1 1],A,b,[],...
    optimset('TolFun',1e-6,'MaxFunEvals',1000000,'TolX',1e-6));