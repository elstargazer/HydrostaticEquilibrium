ccc

a=3;
b=2;
c=1;

rho=1;

<<<<<<< HEAD
step = 0.0005;
x=(0:step:3.5)+0.001;
z=(0:step:3.5)+0.001;
=======
step = 0.0001;
x=3.401:step:3.501;
z=3.401:step:3.501;
>>>>>>> 3a0ebf8cff0cd8c2f9067f3e09b40bc6954ce854

[x,z] = meshgrid(x,z);
y=0.001*z;

angle=linspace(0,2*pi,100);
zell=c*sin(angle);
xell=a*cos(angle);


ax=zeros(size(x));
ay=zeros(size(x));
az=zeros(size(x));
U=zeros(size(x));

tic
parfor i=1:numel(x); 
%     U(i)=GetK(a,b,c,x(i),y(i),z(i))
%     [ax(i),ay(i),az(i)]=GetDK(a,b,c,x(i),y(i),z(i))
    
    [ax(i),ay(i),az(i)] = Ell3Acc(a,b,c,x(i),y(i),z(i),rho);
    U(i) = Ell3Pot(a,b,c,x(i),y(i),z(i),rho);
end
toc

% figure; hold on;
% quiver(x,z,ax,az);
% axis equal

figure; hold on;
<<<<<<< HEAD
quiver(x,z,ax,az,2,'Color','r');
contour(x,z,U);
plot(xell,zell,'-k');
xlim([-3.5 3.5]);
ylim([-3.5 3.5]);
box on;
=======
contour(x,z,U);
plot(xell,zell,'-k');
>>>>>>> 3a0ebf8cff0cd8c2f9067f3e09b40bc6954ce854
axis equal

[fx,fz]=gradient(U,step,step);

<<<<<<< HEAD

=======
>>>>>>> 3a0ebf8cff0cd8c2f9067f3e09b40bc6954ce854
figure; hold on;
quiver(x,z,fx,fz,'Color','b');
quiver(x,z,ax,az,'Color','r');
plot(xell,zell,'-k');
axis equal
<<<<<<< HEAD
xlim([-3.5 3.5]);
ylim([-3.5 3.5]);
=======
xlim([3.4 3.5]);
ylim([3.4 3.5]);
>>>>>>> 3a0ebf8cff0cd8c2f9067f3e09b40bc6954ce854

dv=sqrt((ax-fx).^2+(az-fz).^2);
amag=sqrt(ax.*ax+az.*az);
fmag=sqrt(fx.*fx+fz.*fz);

phi = acos((ax.*fx + az.*fz)./...
    sqrt((ax.*ax+az.*az).*(fx.*fx+fz.*fz)))*180/pi;

figure; hold on;
surf(x,z,log10(abs(dv./fmag))); shading flat
plot(xell,zell,'-k');
axis equal
caxis([-16 5]);
colorbar
<<<<<<< HEAD
xlim([0 3.5]);
ylim([0 3.5]);
=======
xlim([3.4 3.5]);
ylim([3.4 3.5]);
>>>>>>> 3a0ebf8cff0cd8c2f9067f3e09b40bc6954ce854

w = x.*x./a./a + z.*z./c./c;

figure; hold on;
hist(log10(abs(dv(w>1)./fmag(w>1))),30);
box on;
xlim([-8 2]);

% figure; hold on;
% pcolor(x,z,real(phi)); shading flat
% plot(xell,zell,'-k');
% axis equal
% colorbar
% xlim([0 3.5]);
% ylim([0 3.5]);

