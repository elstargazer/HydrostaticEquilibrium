function [J2,C22]=rf2j2c22(r,fq,fp,Rref)

% rrat=r.*r./(Rref.*Rref);
% 
% C22=0.05*(fq.*(fq-2))./((((fq-1).*(fp-1)).^(2/3))).*rrat;
% J2=-0.1*(2*fp.*(fp-2)-fq.*(fq-2))./(((fq-1).*(fp-1)).^(2/3)).*rrat;
% 

[a,b,c]=fr2abc(r,fp,fq);

J2  =-0.20*(c.*c-0.5*(a.*a+b.*b))/(Rref.^2)*NormCoef(2,0);
C22 = 0.05*(a.*a-b.*b)/(Rref.^2)*NormCoef(2,2);

