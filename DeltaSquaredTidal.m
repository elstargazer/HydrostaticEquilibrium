function d2=DeltaSquaredTidal(f,r1,T,rho1)

W=2*pi/(3600*T);

fp1=f(1);
fq1=f(2);

[a1,b1,c1]=fr2abc(r1,fp1,fq1);


% [rho1 rho2]'
% [a1 b1 c1]'
% [a2 b2 c2]'
% [fq2 fp2]'

% 1 outer
% 2 inner

d2=(EllPotTot3Ax(a1,b1,c1,0,b1,0,rho1,W)-...
    EllPotTot3Ax(a1,b1,c1,a1,0,0,rho1,W)).^2+...
    (EllPotTot3Ax(a1,b1,c1,0,0,c1,rho1,W)-...
    EllPotTot3Ax(a1,b1,c1,a1,0,0,rho1,W)).^2;