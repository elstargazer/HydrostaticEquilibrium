% ccc

%% Input parameters
G=6.67e-11;

r1=261110;
r2=145000;
rho1=3157.68;
rho2=5000;

T0=2.8;
T1=17.2;
Tstep=0.01;


% period
T=T1:-Tstep:T0;
omega=2*pi./(T.*3600);

M1=4/3*pi*rho1*r1.^3;
M2=4/3*pi*(rho2-rho1)*r2.^3;

M=M1+M2;

mu=G*M


%% Initial guess

f1i=0.0:0.005:0.9;
f2i=0.0:0.005:0.9;

[f1i,f2i]=meshgrid(f1i,f2i);

progressbar(0);

clear d2

for i=1:numel(f1i)    
    d2(i)=DeltaSquared2l([f1i(i) f2i(i)],r1,r2,T(1),rho1,rho2);
    progressbar(i/numel(f1i));
end

progressbar(1);

d2=reshape(d2,size(f1i));


figure; hold on;

pcolor(f1i,f2i,d2);
contour(f1i,f2i,log10(d2),60,'Color','k');
shading interp


% Get coordinated of the minimum

mincoord=ginput(1);

% Find initial equilibrium

[f0,fval]=HydrostaticStateExact2l(r1,r2,T(1),rho1,rho2,mincoord(1),mincoord(2));

plot(f0(1),f0(2),'*w');

ylim([0 0.9])
xlim([0 0.9])



set(gca,'FontSize',20);

cbar=colorbar;

set(cbar,'FontSize',20);

ylabel(cbar,'log_{10}(\Delta^{2}) [m^{4} sec ^{-4}]','FontSize',20)

xlabel('Outer shape flattening []','FontSize',20);
ylabel('Core flattening []','FontSize',20);

title(['T = ' num2str(T1) ' [hr]'],'FontSize',20);

box on;
%% Main loop

clear d2

progressbar(0);

tic
for i=1:numel(T)
       
    [fh(i,:),d2(i)]=HydrostaticStateExact2l(r1,r2,T(i),rho1,rho2,f0(1),f0(2));
    f0=fh(i,:);
    
    progressbar(i/numel(T));
    
end
toc

progressbar(1);


%% Plot
try
figure(figuref)
catch
end


figure
hold on;

plot(T,fh(:,1),'r-','LineWidth',2);
plot(T,fh(:,2),'b-','LineWidth',2);

NorthFlattening=0.1471;
GlobalFlattening=0.1931;

line([2.8 10.2],[NorthFlattening NorthFlattening],'Color','m');
line([2.8 10.2],[GlobalFlattening GlobalFlattening],'Color','k');

xlim([2.8 10.2])
ylim([0 0.5])

xlabel('Rotation period [hr]','FontSize',12);
ylabel('Flattening = (a-c)/a','FontSize',12);

set(gcf, 'Units','centimeters', 'Position',[0 0 13 9])
set(gcf, 'PaperPositionMode','auto')


%% Checking Radau-Darwin relation
% 1
[a1,c1]=f2axes(r1,fh(:,1));
[a2,c2]=f2axes(r2,fh(:,2));

Ch1=0.2*(M1)*(a1.^2+c1.^2);
Ch2=0.2*(M2)*(a2.^2+c2.^2);

Ch=Ch1+Ch2;

lambdah=Ch./(M.*r1.^2);

q=(omega.^2).*(a1'.^3)./(G.*M);
epsilon=fh(:,1);
eta=5*q'./(2.*epsilon)-2;
lambda_RD=2/3*(1-2/5*sqrt(1+eta));

figure
hold on;
set(gca,'FontSize',20)

plot(T,lambdah,'r-','LineWidth',2);
plot(T,lambda_RD,'b-','LineWidth',2);

xlabel('Rotation period [hr]','FontSize',20);
ylabel('\lambda []','FontSize',20);

% 2

f1_RD=5/2*q./(1+25/4*(1-1.5*lambdah').^2);

figure
hold on;
set(gca,'FontSize',20)

plot(T,fh(:,1),'r-','LineWidth',2);
plot(T,f1_RD,'b-','LineWidth',2);

xlabel('Rotation period [hr]','FontSize',20);
ylabel('Flattening = (a-c)/a','FontSize',20);


% figure; hold on
% data=load('SixOrderHydroList');
% 
% plot(data(3,:),data(1,:),'g');
% plot(data(3,:),data(2,:),'r');


% figure;
% semilogy(T,d2,'b');
% grid on;
% 





