function d2=DeltaSquared2lCoreNonEq(f,r1,r2,T,rho1,rho2,fp1,fq1)

W=2*pi/(3600*T);

fp2=f(1);
fq2=f(2);

[a1,b1,c1]=fr2abc(r1,fp1,fq1);
[a2,b2,c2]=fr2abc(r2,fp2,fq2);

% 1 outer
% 2 inner

<<<<<<< HEAD
% d2=(EllPotTot2l(a1,c1,a2,c2,0,0,c1,rho1,rho2,W)-EllPotTot2l(a1,c1,a2,c2,a1,0,0,rho1,rho2,W)).^2+...
%   (EllPotTot2l(a1,c1,a2,c2,0,0,c2,rho1,rho2,W)-EllPotTot2l(a1,c1,a2,c2,a2,0,0,rho1,rho2,W)).^2;

% d2=(EllPotTot3Ax2l(a1,b1,c1,a2,b2,c2,0,b1,0,rho1,rho2,W)-...
%     EllPotTot3Ax2l(a1,b1,c1,a2,b2,c2,a1,0,0,rho1,rho2,W)).^2+...
%    (EllPotTot3Ax2l(a1,b1,c1,a2,b2,c2,0,0,c1,rho1,rho2,W)-...
%    EllPotTot3Ax2l(a1,b1,c1,a2,b2,c2,a1,0,0,rho1,rho2,W)).^2+...
%    (EllPotTot3Ax2l(a1,b1,c1,a2,b2,c2,0,b2,0,rho1,rho2,W)-...
%    EllPotTot3Ax2l(a1,b1,c1,a2,b2,c2,a2,0,0,rho1,rho2,W)).^2+...
%    (EllPotTot3Ax2l(a1,b1,c1,a2,b2,c2,0,0,c2,rho1,rho2,W)-...
%    EllPotTot3Ax2l(a1,b1,c1,a2,b2,c2,a2,0,0,rho1,rho2,W)).^2;

try
d2=(EllPotTot3Ax2l(a1,b1,c1,a2,b2,c2,0,b1,0,rho1,rho2,W)-...
    EllPotTot3Ax2l(a1,b1,c1,a2,b2,c2,a1,0,0,rho1,rho2,W)).^2+...
    (EllPotTot3Ax2l(a1,b1,c1,a2,b2,c2,0,0,c1,rho1,rho2,W)-...
    EllPotTot3Ax2l(a1,b1,c1,a2,b2,c2,a1,0,0,rho1,rho2,W)).^2+...
    (EllPotTot3Ax2l(a1,b1,c1,a2,b2,c2,0,b1,0,rho1,rho2,W)-...
    EllPotTot3Ax2l(a1,b1,c1,a2,b2,c2,0,0,c1,rho1,rho2,W)).^2;
catch
    
    1
    
end
% 
=======
d2=(EllPotTot3Ax2l(a1,b1,c1,a2,b2,c2,0,b1,0,rho1,rho2,W)-...
    EllPotTot3Ax2l(a1,b1,c1,a2,b2,c2,a1,0,0,rho1,rho2,W)).^2+...
   (EllPotTot3Ax2l(a1,b1,c1,a2,b2,c2,0,0,c1,rho1,rho2,W)-...
   EllPotTot3Ax2l(a1,b1,c1,a2,b2,c2,a1,0,0,rho1,rho2,W)).^2;

>>>>>>> 3a0ebf8cff0cd8c2f9067f3e09b40bc6954ce854
%% Miminize globally
% 
% TRr=IcosahedronMesh;
% TR_2r=SubdivideSphericalMesh(TRr,2);
% 
% [lambda,fi,~]=cart2sph(TR_2r.X(:,1),TR_2r.X(:,2),TR_2r.X(:,3));
% 
% % AGUaxes
% % plotm(fi*180/pi,lambda*180/pi,'.');
% 
% [x1,y1,z1]=TriEllRadVec(fi,lambda,a1,b1,c1,'xyz');
% % [x2,y2,z2]=TriEllRadVec(fi,lambda,a2,b2,c2,'xyz');
% 
% U1=zeros(size(x1));
% % U2=U1;
% 
% for i=1:numel(x1)  
%     U1(i)=EllPotTot3Ax2l(a1,b1,c1,a2,b2,c2,x1(i),y1(i),z1(i),rho1,rho2,W);  
% end
% 
% % for i=1:numel(x2)
% %     U2(i)=EllPotTot3Ax2l(a1,b1,c1,a2,b2,c2,x2(i),y2(i),z2(i),rho1,rho2,W);
% % end
% % 
% % d2=range(U1).^2 + range(U2).^2;
% 
% d2=range(U1).^2;
% 
% AGUaxes
% scatterm(fi*180/pi,lambda*180/pi,100,U1,'filled');


