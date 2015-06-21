function [J2,C22]=RadFlat2J2Tri(r1,r2,fp1,fp2,fq1,fq2,rho1,rho2,Rref)

% [a1,b1,c1]=fr2abc(r1,fp1,fq1);
% [a2,b2,c2]=fr2abc(r2,fp2,fq2);

V1=4/3*pi.*r1.^3;
V2=4/3*pi.*r2.^3;

M1=rho1.*V1;
M2=(rho2-rho1).*V2;

M=M1+M2;
% 
% j2 = @(a,b,c, Rref) -0.2*(c.*c-0.5*(a.*a+b.*b))./(Rref.*Rref);
% c22 = @(a,b,Rref) 0.05*(a.*a-b.*b)./(Rref.*Rref);
% 
% J2=j2(a1,b1,c1,Rref).*M1+j2(a2,b2,c2,Rref).*M2;
% J2=J2./M;
% 
% C22=c22(a1,b1,Rref).*M1+c22(a2,b2,Rref).*M2;
% C22=C22./M;


[J2_1,C22_1]=rf2j2c22(r1,fq1,fp1,Rref);
[J2_2,C22_2]=rf2j2c22(r2,fq2,fp2,Rref);

J2=(J2_1.*M1+J2_2.*M2)./M;
C22=(C22_1.*M1+C22_2.*M2)./M;