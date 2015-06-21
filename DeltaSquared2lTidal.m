function d2=DeltaSquared2lTidal(f,r1,r2,T,rho1,rho2)

W=2*pi/(3600*T);

fp1=f(1);
fp2=f(2);
fq1=f(3);
fq2=f(4);

[a1,b1,c1]=fr2abc(r1,fp1,fq1);
[a2,b2,c2]=fr2abc(r2,fp2,fq2);

% [rho1 rho2]'
% [a1 b1 c1]'
% [a2 b2 c2]'
% [fq2 fp2]'

% 1 outer
% 2 inner

d2=(EllPotTot3Ax2l(a1,b1,c1,a2,b2,c2,0,b1,0,rho1,rho2,W)-...
    EllPotTot3Ax2l(a1,b1,c1,a2,b2,c2,a1,0,0,rho1,rho2,W)).^2+...
   (EllPotTot3Ax2l(a1,b1,c1,a2,b2,c2,0,0,c1,rho1,rho2,W)-...
   EllPotTot3Ax2l(a1,b1,c1,a2,b2,c2,a1,0,0,rho1,rho2,W)).^2+...
   (EllPotTot3Ax2l(a1,b1,c1,a2,b2,c2,0,b2,0,rho1,rho2,W)-...
   EllPotTot3Ax2l(a1,b1,c1,a2,b2,c2,a2,0,0,rho1,rho2,W)).^2+...
   (EllPotTot3Ax2l(a1,b1,c1,a2,b2,c2,0,0,c2,rho1,rho2,W)-...
   EllPotTot3Ax2l(a1,b1,c1,a2,b2,c2,a2,0,0,rho1,rho2,W)).^2;