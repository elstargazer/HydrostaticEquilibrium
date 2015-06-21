function d2=DeltaSquared3Ax2lmod(x,r1,r2,T,fq1,fp1)

W=2*pi/(3600*T);

rho1=x(1);
rho2=x(2);
fq2=x(3);
fp2=x(4);


a1=r1/(((fp1-1)*(fq1-1))^(1/3));
b1=(r1 - fq1*r1)/(((fp1-1)*(fq1-1))^(1/3));
c1=(r1 - fp1*r1)/(((fp1-1)*(fq1-1))^(1/3));

a2=r2/(((fp2-1)*(fq2-1))^(1/3));
b2=(r2 - fq2*r2)/(((fp2-1)*(fq2-1))^(1/3));
c2=(r2 - fp2*r2)/(((fp2-1)*(fq2-1))^(1/3));

% [rho1 rho2]'
% 
% [a1 b1 c1]'
% [a2 b2 c2]'
% 
% [fq2 fp2]'


% 1 outer
% 2 inner

d2=(EllPotTot3Ax2l(a1,b1,c1,a2,b2,c2,0,b1,0,rho1,rho2,W)-EllPotTot3Ax2l(a1,b1,c1,a2,b2,c2,a1,0,0,rho1,rho2,W)).^2+...
   (EllPotTot3Ax2l(a1,b1,c1,a2,b2,c2,0,0,c1,rho1,rho2,W)-EllPotTot3Ax2l(a1,b1,c1,a2,b2,c2,a1,0,0,rho1,rho2,W)).^2+...
   (EllPotTot3Ax2l(a1,b1,c1,a2,b2,c2,0,b2,0,rho1,rho2,W)-EllPotTot3Ax2l(a1,b1,c1,a2,b2,c2,a2,0,0,rho1,rho2,W)).^2+...
   (EllPotTot3Ax2l(a1,b1,c1,a2,b2,c2,0,0,c2,rho1,rho2,W)-EllPotTot3Ax2l(a1,b1,c1,a2,b2,c2,a2,0,0,rho1,rho2,W)).^2;