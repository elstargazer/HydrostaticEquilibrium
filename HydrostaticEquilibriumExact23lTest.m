%% Input parameters


r1=261100;
r2=110000;
rho1=3107;
rho2=7800;

T0=2;
T1=10;
Tstep=0.01;


% period
T=T1:-Tstep:T0;


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

plot(T,fh(:,1),'r-.','LineWidth',2);
plot(T,fh(:,2),'b-.','LineWidth',2);

% figure; hold on
% data=load('SixOrderHydroList');
% 
% plot(data(3,:),data(1,:),'g');
% plot(data(3,:),data(2,:),'r');


figure;
semilogy(T,d2,'b');
grid on;
% 





