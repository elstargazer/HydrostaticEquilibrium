function [J2,C22]=RadFlat2J2C22(r1,r2,fp1,fq1,fp2,fq2,rho1,rho2,Rref)


[a1,b1,c1]=fr2abc(r1,fp1,fq1);
[a2,b2,c2]=fr2abc(r2,fp2,fq2);


V1=4/3*pi.*r1.^3;
V2=4/3*pi.*r2.^3;

M1=rho1.*V1;
M2=(rho2-rho1).*V2;

M=M1+M2;

J2 = -((c1.*c1-0.5*(a1.*a1+b1.*b1)).*M1+...
    (c2.*c2-0.5*(a2.*a2+b2.*b2)).*M2)./...
    (5*Rref*Rref.*M)*NormCoef(2,0);

C22 = ((a1.*a1-b1.*b1).*M1 + ...
    (a2.*a2-b2.*b2).*M2)./...
    (20*Rref*Rref.*M)*NormCoef(2,2);
