function PlotPot2El(r1,r2,fp1,fq1,fp2,fq2,rho1,rho2,W)

fi=(-90:2:90)/180*pi;
lambda=(-180*2:180)/180*pi;

[fi,lambda] = meshgrid(fi, lambda);

[a1,b1,c1]=fr2abc(r1,fp1,fq1);
[a2,b2,c2]=fr2abc(r2,fp2,fq2);

[x1,y1,z1]=TriEllRadVec(fi,lambda,a1,b1,c1,'xyz');

U = zeros(size(fi));

progressbar(0);

parfor i=1:numel(fi)
    U(i) = EllPotTot3Ax2l(a1,b1,c1,a2,b2,c2,x1(i),y1(i),z1(i),rho1,rho2,W);
%     progressbar(i/numel(fi));
end

progressbar(1);

U = reshape(U,size(fi));

U_mean = mean(U(:));

g = 0.28;
dr = (U - U_mean)/g;

AGUaxes;
pcolorm(fi*180/pi,lambda*180/pi,dr); shading interp;

cbar = colorbar;
ylabel(cbar,'Geoid Height Anomaly [m]','FontSize',20);