% ccc
% 
% M=1;
% R=1;
% 
% G=6.67e-11;
% 
% x=-20:0.05:20;
% z=-20:0.05:20;
% 
% [x,z]=meshgrid(x,z);
% 
% y=0;
% 
% r=sqrt(x.*x+y.*y+z.*z);
% 
% U=RingPot(x,y,z,R,M);
% Ur=G.*M./r;
% 
% 
% figure;
% hold on;
% 
% pcolor(x,z,U); shading interp;
% axis equal
% colorbar
% 
% figure;
% hold on;
% 
% pcolor(x,z,Ur); shading interp;
% axis equal
% colorbar
% 
% figure;
% hold on;
% 
% plot(r(200,:),log10(abs(U(200,:))),'-r');
% plot(r(200,:),log10(abs(Ur(200,:))),'-b');
% 
% figure;
% hold on;
% 
% plot(r(200,:),log10(abs((Ur(200,:)-U(200,:))./Ur(200,:))),'-b');

ccc

r=6371000;
rho=5517;


T0=4;
T1=7*24;
Tstep=0.5;

% Rr=17536000;

Rr=r*3;

% Mr=1.52e22; % Charon mass
Mr=7.34e22; % Moon mass


T=T1:-Tstep:T0;

f0=0;

fh=zeros(numel(T),1);
d2=zeros(numel(T),1);

tic
for i=1:numel(T)       
    [fh(i),d2(i)]=HydrostaticStateExact(r,T(i),rho,0);
    f0=fh(i);    
end
toc

fhr=zeros(numel(T),1);
d2r=zeros(numel(T),1);


progressbar(0);
tic
parfor i=1:numel(T)       
    [fhr(i),d2r(i)]=HydrostaticStateExactRing(r,T(i),rho,Rr,Mr,0.1);
    f0=fhr(i);
%     progressbar(i/numel(T));
end
toc
progressbar(1);



figure; hold on;
plot(T,fhr,'g-');
plot(T,fh,'b-');


% figure; hold on;
% plot(T,fhr-fh,'b-');




%%

ccc;


M=4.006e21;

% T=3.915341;
T=22

a_guess=960000;
b_guess=770000;
c_guess=495000;

r1=(a_guess*b_guess*c_guess)^(1/3);
V=4/3*pi*r1^3;

rhomean=M/V;

fstep=0.005;

maxfq=0.99;
maxfp=0.99;

fqi=0.00001:fstep:maxfq;
fpi=0.00001:fstep:maxfp;

[fqi,fpi]=meshgrid(fpi,fqi);


for i=1:numel(fpi)
      if (fpi(i)>fqi(i))
         d2(i)=DeltaSquared3AxRing([fpi(i) fqi(i)],r1,T,rhomean);       
      else
        d2(i)=NaN;
      end    
%          progressbar(i/numel(fpi));
end

d2=reshape(d2,size(fpi));

figure; hold on;
set(gca,'FontSize',20);


pcolor(fpi,fqi,log10(d2));
contour(fpi,fqi,log10(d2),60,'Color','k'); shading interp;

cbar=colorbar;
ylabel(cbar,'log_{10}(\Delta^{2}) [m^{4} sec ^{-4}]','FontSize',20)

xlim([0 1]);
ylim([0 1]);

xlabel('Polar flattening []','FontSize',20); 
ylabel('Equatorial flattening []','FontSize',20); 
caxis([6 16]);







