function [fh,fval]=HydrostaticStateExact3Ax(r,T,rho,fp0,fq0)

A=[-1 1];

b=[0];

[fh,fval]=fminsearchcon(@(f) DeltaSquared3Ax(f,r,T,rho),...
    [fp0 fq0],[0 0],[1 1],A,b,[],...
    optimset('TolFun',1e-10,'MaxFunEvals',10000));



