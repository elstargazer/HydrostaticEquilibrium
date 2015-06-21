function [fh,fval]=HydrostaticStateExact3l(r1,r2,r3,T,rho1,rho2,rho3,f10,f20,f30)


[fh,fval]=fminsearch(@(f) DeltaSquared3l(f,r1,r2,r3,T,rho1,rho2,rho3),[f10 f20 f30],...
    optimset('TolFun',1e-20));




