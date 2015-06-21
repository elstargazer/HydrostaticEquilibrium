 ccc

%% Input parameters

GM=17.28e9;
G=6.67e-11;
M=GM/G;

V=0.74886e8*1e9;

h=40000; % crustal thickness


r1=(V/(4/3*pi))^(1/3);
r2=r1-h;
r3=114200;

rho3=6259;
rho1=2885;



% ro_crust=2885;
% ro_mantle=3266;
% ro_core=6259;




V=4/3*pi*r1^3;

V1=4/3*pi*r1^3;
V2=4/3*pi*r2^3;
V3=4/3*pi*r3^3;

M3t=V3*(rho3);
M1t=(V1-V2)*(rho1);
M2t=M-M1t-M3t;

V2t=V2-V3;

rho2=M2t/V2t


T0=2.8;
T1=10.2;
Tstep=0.01;


% period
T=T1:-Tstep:T0;
%%

M1=V1*rho1;
M2=V2*(rho2-rho1);
M3=V3*(rho3-rho2);

Mc=M1+M2+M3;


%% Initial guess

% f1i=0.0:0.005:0.9;
% f2i=0.0:0.005:0.9;
% 
% 
% [f1i,f2i]=meshgrid(f1i,f2i);
% 
% progressbar(0);
% 
% clear d2
% 
% for i=1:numel(f1i)    
%     d2(i)=DeltaSquared2l([f1i(i) f2i(i)],r1,r2,T(1),rho1,rho2);
%     progressbar(i/numel(f1i));
% end
% 
% progressbar(1);
% 
% d2=reshape(d2,size(f1i));
% 
% 
% figure; hold on;
% 
% pcolor(f1i,f2i,d2);
% contour(f1i,f2i,log10(d2),60,'Color','k');
% shading interp
% 
% 
% % Get coordinated of the minimum
% 
% mincoord=ginput(1);

% Find initial equilibrium

f10=0.1;
f20=0.1;
f30=0.1;

[f0,fval]=HydrostaticStateExact3l(r1,r2,r3,T(1),rho1,rho2,rho3,f10,f20,f30);



%% Main loop

clear d2

progressbar(0);

tic
for i=1:numel(T)
       
    [fh(i,:),d2(i)]=HydrostaticStateExact3l(r1,r2,r3,T(i),rho1,rho2,rho3,f0(1),f0(2),f0(3));
    f0=fh(i,:);
    
    progressbar(i/numel(T));
    
end
toc

progressbar(1);


%% Plot
% try
% figure(figuref)
% catch
% end

figurevf=figure;
set(gca,'FontSize',12);

hold on;

plot(T,fh(:,1),'r-','LineWidth',2);
plot(T,fh(:,2),'-','LineWidth',1,'Color',[0 0.6 0]);
plot(T,fh(:,3),'b-','LineWidth',2);

CurrentRotationPeriod=5.342;
HydrostaticRotationPeriod=4.831;
line([CurrentRotationPeriod CurrentRotationPeriod],[0 1],'Color','k')
line([HydrostaticRotationPeriod HydrostaticRotationPeriod],[0 1],'Color','k')

line([2.77 2.77],[0 1],'Color','b')


NorthFlattening=0.1471;
GlobalFlattening=0.1931;

line([2.0 10.2],[NorthFlattening NorthFlattening],'Color','k');
line([2.0 10.2],[GlobalFlattening GlobalFlattening],'Color','k');


xlim([2.8 10.2])
ylim([0 0.6])

xlabel('Rotation period [hr]','FontSize',12);
ylabel('Flattening = (a-c)/a','FontSize',12);



set(gcf, 'Units','centimeters', 'Position',[0 0 13 9])
set(gcf, 'PaperPositionMode','auto')

legend({'Outer shape flattening','Mantle flattening','Core flattening'});


set(gcf, 'Units','centimeters', 'Position',[0 0 13 9])
set(gcf, 'PaperPositionMode','auto')

box on


print(figurevf, '-dpsc', 'HydrostaticEquilibriumFigureExact.eps');


% figure; hold on
% data=load('SixOrderHydroList');
% 
% plot(data(3,:),data(1,:),'g');
% plot(data(3,:),data(2,:),'r');


figure; hold on;

set(gca,'FontSize',20);
xlabel('Rotation period [hr]','FontSize',20);
ylabel('log_{10}(\Delta^{2}) [m^{4} sec ^{-4}]','FontSize',20)
box on

plot(T,d2,'b');
set(gca,'YScale','Log')
grid on;
xlim([2.0 10.2])






