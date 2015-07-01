function [fh,fval]=HydrostaticStateExact2l(r1,r2,T,rho1,rho2,f10,f20)

[fh,fval]=fminsearch(@(f) DeltaSquared2l(f,r1,r2,T,rho1,rho2),[f10 f20],...
    optimset('TolFun',1e-10));




