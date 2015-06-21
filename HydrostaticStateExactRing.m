function [fh,fval]=HydrostaticStateExactRing(r,T,rho,Rr,Mr,f0)


[fh,fval]=fminsearch(@(f) DeltaSquaredRing(f,r,T,rho,Rr,Mr),f0,...
    optimset('TolFun',1e-10));

%% plot potential
% W=2*pi/(3600*T);
% 
% a=r/((1-fh).^(1/3));
% c=r/((1-fh).^(1/3))-fh*r/((1-fh).^(1/3));
% 
% fi=(-90:1:90)/180*pi;
% lambda=(-180:1:180)/180*pi;
% 
% [fi,lambda]=meshgrid(fi,lambda);
% 
% [x,y,z]=TriEllRadVec(fi,lambda,a,a,c);
% 
% 
% U=EllPotTotRing(a,c,x,y,z,rho,W,Rr,Mr);

% figure; hold on;
% pcolor(lambda*180/pi,fi*180/pi,U); shading interp;
% colorbar;

