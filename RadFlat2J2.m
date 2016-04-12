function J2=RadFlat2J2(r1,r2,f1,f2,rho1,rho2,Rref)

a1=r1./((1-f1).^(1/3));
c1=r1./((1-f1).^(1/3))-f1.*r1./((1-f1).^(1/3));

a2=r2./((1-f2).^(1/3));
c2=r2./((1-f2).^(1/3))-f2.*r2./((1-f2).^(1/3));

V1=4/3*pi.*r1.^3;
V2=4/3*pi.*r2.^3;

M1=rho1.*V1;
M2=(rho2-rho1).*V2;

M=M1+M2;


J2_1 = - (c1.*c1-a1.*a1)/(5*Rref*Rref)/sqrt(5);
J2_2 = - (c2.*c2-a2.*a2)/(5*Rref*Rref)/sqrt(5);

J2=(J2_1.*M1+J2_2.*M2)./M;

% 
% disp('J2_1');
% -(c1.*c1-a1.*a1)./(5*Rref*Rref)/sqrt(5)
% 
% disp('J2_2');
% -(c2.*c2-a2.*a2)./(5*Rref*Rref)/sqrt(5)

