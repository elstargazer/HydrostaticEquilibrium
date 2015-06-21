function [fh,fval]=HydrostaticStateExact(r,T,rho,f0)


[fh,fval]=fminsearch(@(f) DeltaSquared(f,r,T,rho),f0,...
    optimset('TolFun',1e-10));




