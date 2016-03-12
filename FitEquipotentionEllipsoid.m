% function [J2,C22] = FitEquipotentionEllipsoid(ell)
ccc

T=9.073859324514187; % DLR
% ell = [484.17 481.41 447.80]*1000;
ell = [483.43 481.36 445.68]*1000; % PreSurvey SPC JPL

a=ell(1);
b=ell(2);
c=ell(3);

Rref=476000;
GM=62.68e9;

x=a;
y=b;
z=c;

% J2o=0.0001;
% C22o=0.00001;
% 
% lmcosi = [0 0 1 0;
%           1 0 0 0;
%           1 1 0 0;
%           2 0 -J2o 0;
%           2 1 0 0;
%           2 2 C22o 0];
% 
% U0 = GM/sqrt(x.*x+y.*y+z.*z);
% 
% U1=GravityPotential(GM,Rref,lmcosi,x,y,z);
% U2=SecondDegreePot([J2o C22o],GM,x,y,z,Rref);
% 
% [U1-U0 U2-U0]
% 
% (U1-U2)/(U1)

Delta2 = @(JS) (SecondDegreePot(JS,GM,a,0,0,Rref,T)-...
                SecondDegreePot(JS,GM,0,b,0,Rref,T)).^2+...
               (SecondDegreePot(JS,GM,0,b,0,Rref,T)-...
                SecondDegreePot(JS,GM,0,0,c,Rref,T)).^2+...
               (SecondDegreePot(JS,GM,a,0,0,Rref,T)-...
                SecondDegreePot(JS,GM,0,0,c,Rref,T)).^2;
            
Delta2b = @(JS) (SecondDegreePot(JS,GM,a,0,0,Rref,T)-...
                SecondDegreePot(JS,GM,0,b,0,Rref,T)).^2+...
               (SecondDegreePot(JS,GM,a,0,0,Rref,T)-...
                SecondDegreePot(JS,GM,0,0,c,Rref,T)).^2;
           

[JS,fval]=fminsearch(Delta2, [0 0],...
    optimset('MaxFunEvals',1e4,'TolX',1e-10));

[JSb,fval]=fminsearch(Delta2b, [0 0],...
    optimset('MaxFunEvals',1e4,'TolX',1e-10));

J2=JS(1);
C22=JS(2);

%% Plot potential
% lambda=(-180:1:180)/180*pi;
% fi=(-90:1:90)/180*pi;
% [lambda,fi]=meshgrid(lambda,fi);
% [x,y,z]=TriEllRadVec(fi,lambda,a,b,c,'xyz');
% 
% U=SecondDegreePot(JS,GM,x,y,z,Rref);
% 
% 
% AGUaxes;
% pcolorm(fi*180/pi,lambda*180/pi,U);
% colorbar;
% 
% range(U(:))./mean(U(:))


J2

C22

fval