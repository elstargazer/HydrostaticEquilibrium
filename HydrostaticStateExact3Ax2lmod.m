function [xh,fval,exitflag,output]=HydrostaticStateExact3Ax2lmod(fq1,fp1,r1,r2,T,x_guess)


% [fh,fval]=fminsearch(@(f) DeltaSquared3Ax2l(f,r1,r2,T,rho1,rho2),f_guess,...
%     optimset('TolFun',1e-10));

A=[1 -1 0 0;
   0 0 1 -1];

b=[0 0];

% DeltaSquared3Ax2lmod(x_guess,r1,r2,T,fq1,fp1)


 [xh,fval,exitflag,output]=fminsearchcon(@(x) DeltaSquared3Ax2lmod(x,r1,r2,T,fq1,fp1),...
     x_guess,[500 1000 0 0],[4000 5000 1 1],A,b,[],...
     optimset('TolFun',1e-6,'MaxFunEvals',10000000,'TolX',1e-4));