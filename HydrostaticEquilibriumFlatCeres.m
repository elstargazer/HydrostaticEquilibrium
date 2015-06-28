% close all
ccc;
tic;
%% Input parameters
G=6.67384e-11;

GM=62.68e9; % PCK version 0.2
T=9.073859324514187; % DLR

Rref=500000;

% Equatorial semi-axis (km)   487.30  1.8
% Polar semi-axis (km)        454.70  1.6

% a=[487.3 479.7 483.4990 482.1130 482.7880];
% c=[454.7 444.4 447.6 446.1255 447.8000];
% sigmaa=[1.8 2.3 0.7 0 0];
% sigmac=[1.6 2.1 0.4 0 0];

a=[482.1130 482.7880];

% OpNav5 JPL shape
a3=484.17;
b3=481.41;
c3=447.80;

fpa_obs=(a3-c3)/a3;
fq_obs=(a3-b3)/a3;
fpb_obs=(b3-c3)/b3;

fp_obs=(sqrt(a3*b3)-c3)/sqrt(a3*b3);

c=[446.1255 447.8000];
sigmaa=[0.0 0.0];
sigmac=[0.0 0.0];

Nconf=numel(a);
ModelIn=1;

%% figure settings 

fntsize = 12;
fntsize_sm = 10;
im_size=[0 0 13 9];

fig_folder='~/Dawn/Papers/CeresPaper1/';

%% Preliminary conmutations

M=GM/G;

R=(a.*a.*c).^(1/3);

V=1e9*4/3*pi*R.^3;

fouter_obs=(a-c)./a;
sigma_fouter_obs=sqrt(...
    (c.*c.*sigmaa.*sigmaa+a.*a.*sigmac.*sigmac)./(c.^4));

router=(V./(4/3*pi)).^(1/3);

fcore0=0.1;
fouter0=0.1;

rhomean=M./V;

% router=router(ModelIn);
rhomean=rhomean(ModelIn);
g=GM./(router.^2);
% fouter_obs=fouter_obs(ModelIn);

%% Grid of core radii and densities
% rcore=2500:500:470000;
% rhocoreg=2000:100:4000;

rcore=linspace(10000,470000,100);
rhocoreg=linspace(rhomean,5000,100);

[rhocorei,rcorei]=meshgrid(rhocoreg,rcore);

for i=1:Nconf   
    rhoouteri{i}=-(3*M-4*pi*(rcorei.^3).*rhocorei)./(4*pi*(rcorei.^3)-4*pi*(router(i)^3));
    rhoouteri{i}(rhoouteri{i}<0)=NaN;  
    [fcorei{i},fouteri{i}]=HydrostaticStateExact2lGrid(router(i),rcorei,T,rhoouteri{i},rhocorei);
end

%% Plotting settings

figure;
set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
hold on;box on;grid on;

xlim([10 470]);
ylim([rhomean 5000]);

xlabel('Core size [km]','FontSize',fntsize);
ylabel('Core density [kg/m^3]','FontSize',fntsize);

%% Contour mantle density

pcolor(rcorei/1000,rhocorei,rhoouteri{ModelIn}); shading interp;
cbar=colorbar('FontSize',fntsize);
ylabel(cbar,'Outer density [kg/m^3] ','FontSize',fntsize);

caxis([920 rhomean]);
plot(get(gca,'xlim'),[rhomean rhomean],'--k','LineWidth',3);

%% Contour core flattening

% levels=0:0.01:0.6;
% [C,h]=contour(rcorei/1000,rhocorei,fcorei,levels,'Color',[0 1 0]);
% text_handle = clabel(C,h);

%% Contour outerflattening

fouteri_n=fouteri{ModelIn};
fouteri_n((fouteri_n<0) | isnan(rhoouteri{ModelIn}) )= NaN;

levels=0.045:0.005:0.6;
[C,h]=contour(rcorei/1000,rhocorei,fouteri_n,levels,'Color',...
    [0.7 0.0 0.45],'LineWidth',2);

clabel(C,h,'manual','FontSize',fntsize_sm,'Color','k','EdgeColor','k','BackGroundColor','w');

ccmap=[1 0 0; 0 0 1];

for i=1:Nconf
    
    fouteri_n=fouteri{ModelIn};
    fouteri_n((rhoouteri{ModelIn}<0) | isnan(rhoouteri{ModelIn}) )=NaN;
    
    levels=[fouter_obs(i)-sigma_fouter_obs(i) fouter_obs(i)+sigma_fouter_obs(i)];
    [Chyd{i},h]=contour(rcorei/1000,rhocorei,fouteri_n,levels,'Color',ccmap(i,:),'LineWidth',5,...
        'LineStyle','-');  
end

level_set1=[920 920];
level_set2=[1000 1250 1500 1750 2000];

[C,h]=contour(rcorei/1000,rhocorei,rhoouteri{ModelIn},level_set1,'--',...
    'Color','k','LineWidth',1); 

clabel(C,h,'manual','FontSize',fntsize_sm,...
    'Color','k','EdgeColor','k','BackGroundColor','w');
 
% [C,h]=contour(rcorei/1000,rhocorei,rhoouteri{ModelIn},level_set2,...
%     'Color','c','LineWidth',4);
% clabel(C,h,'manual','FontSize',fntsize_sm,'Color','k','EdgeColor','k','BackGroundColor','w');

cbar=colorbar('FontSize',fntsize);
ylabel(cbar,'Outer density [kg/m^3] ','FontSize',fntsize);


PrintWhite([fig_folder 'Fig_2layer.jpg']);

%% Plot ice shell thickness
figure;
set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
hold on;box on;grid on;

for i=1:Nconf
    
    rcore_h=Chyd{i}(1,:);
    rhocore_h=Chyd{i}(2,:);
    
    rhoouter_h=griddata(rcorei/1000,rhocorei,...
        rhoouteri{i},rcore_h,rhocore_h,'linear');
    plot(rhoouter_h,(router(i)/1000-rcore_h),...
        '-','LineWidth',3,'Color',ccmap(i,:));  
end

xlabel('Shell density [kg/m^{3}]','FontSize',fntsize);
ylabel('Shell thickness [km]','FontSize',fntsize);
xlim([800 rhomean])

legend({'SPG','SPC'},'FontSize',fntsize_sm);
PrintWhite([fig_folder 'Fig_IceThickness.jpg']);


%% Non hydrostatic core

% gi = ginput(1);
% 
% rcore_test = gi(1);
% rho_core_test = gi(2);
% rhoouter_test=-(3*M-4*pi*(rcore_test.^3).*rho_core_test)./...
%     (4*pi*(rcore_test.^3)-4*pi*(router(2)^3));

% [fp_noneq_test,fq_noneq_test,fval,exitflag]=...
%     HydrostaticStateExact2lGridCoreNonEq(...
%     router(2),rcore_test*1000,T,rhoouter_test,...
%     rho_core_test,fpa_obs,fq_obs);

tic
[fp_noneq,fq_noneq,fval,exitflag]=...
    HydrostaticStateExact2lGridCoreNonEq(...
    router(2),rcorei,T,rhoouteri{2},rhocorei,fpa_obs,fq_obs);
toc

% convergence conditions
conv_cond = log10(fval) < 0;

fq_noneq(~conv_cond)=NaN;
fp_noneq(~conv_cond)=NaN;

figure;
set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
hold on;box on;grid on;

xlim([10 450]);
ylim([rhomean 5000]);
caxis([920 rhomean]);

f_ratio = fp_noneq./fq_noneq;

f_levels=0:0.1:1;

pcolor(rcorei/1000,rhocorei,rhoouteri{ModelIn}); shading interp;

[C1,h1] = contour(rcorei/1000,rhocorei,fp_noneq,f_levels,...
    '-','Color','r','ShowText','on');
h1_ = plot(NaN, '-r');

[C2,h2] = contour(rcorei/1000,rhocorei,fq_noneq,f_levels,...
    '-','Color','b','ShowText','on');
h2_ = plot(NaN, '-b');

contour(rcorei/1000,rhocorei,f_ratio,[1 1],'-k',...
    'LineWidth',3,'HandleVisibility','off');
h3_ = plot(NaN, '-k','LineWidth',3);

legend([h1_ h2_,h3_],...
    {'$\Delta f_{p,core}$','$\Delta f_{q,core}$','$f_{q,core}=f_{p,core}$'},...
    'FontSize',fntsize_sm,'Interpreter','latex');

cbar=colorbar('FontSize',fntsize);
ylabel(cbar,'Outer density [kg/m^3] ','FontSize',fntsize);

xlabel('Core size [km]','FontSize',fntsize);
ylabel('Core density [kg/m^3]','FontSize',fntsize);

PrintWhite([fig_folder 'Fig_NonhydroCore.jpg']);

%% Computing J2 and C22 for non-hydrostatic core case

Mcore = 4/3*pi.*(rcorei.^3).*(rhocorei-rhoouteri{2});
Mouter = 4/3*pi.*(router(2).^3).*(rhoouteri{2});

[J2outer,C22outer]=rf2j2c22(router(2),fq_obs,fp_obs,Rref);
[J2core,C22core]=rf2j2c22(rcorei,fq_noneq,fp_noneq,Rref);

J2 = (J2core.*Mcore+J2outer.*Mouter)./M;
C22 = (C22core.*Mcore+C22outer.*Mouter)./M;

figure;
set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
hold on;box on;grid on;

xlim([10 450]);
ylim([rhomean 5000]);
caxis([920 rhomean]);

pcolor(rcorei/1000,rhocorei,rhoouteri{ModelIn}); shading interp;

[C1,h1] = contour(rcorei/1000,rhocorei,J2,20,...
    '-','Color','r','ShowText','on');
h1_ = plot(NaN, '-r');

[C2,h2] = contour(rcorei/1000,rhocorei,C22,20,...
    '-','Color','b','ShowText','on');
h2_ = plot(NaN, '-b');

legend([h1_ h2_],...
    {'$J_{2}$','$C_{22}$'},'FontSize',fntsize_sm,'Interpreter','latex');


PrintWhite([fig_folder 'Fig_GravityPrediction.jpg']);

%% Computing hydrostatic J2
% 
% J2hi=RadFlat2J2(router,rcorei,fouteri,fcorei,rhoouteri,rhocorei,Rref);
% 
% %nmax=10;
% % J=HydroCoeffs2Layer(router,rcorei,fouteri,fcorei,rhoouteri,rhocorei,Rref,nmax);
% 
% 
% levels=0:10e-4:0.04;
% [C,h]=contour(rcorei/1000,rhocorei,J2hi,levels,'Color','m','LineStyle','-','LineWidth',2);
% 
% clabel(C,h,'manual','FontSize',12,'Color','k','EdgeColor','k','BackGroundColor','w');
% 
% % dawn J2 mine
% J2_4=griddata(rcorei/1000,rhocorei,J2hi,C4(1,:),C4(2,:),'cubic');
% J2_4_min=min(J2_4);
% J2_4_max=max(J2_4);
% 
% [J2_4_min J2_4_max]


% [C,h]=contour(rcorei/1000,rhocorei,J2hi,'Color',[0 1 1],'LineStyle','-','LineWidth',2);

% nc=50;
% mc=30;
% 
% for i=1:size(J,3)
%     Jslice(i)=J(nc,mc,i);
% end
% 
% figure; hold on;
% set(gca,'FontSize',20,'YScale','log');
% 
% 
% plot(2:2:nmax,Jslice,'-o','LineWidth',3,'MarkerSize',6);
% 
% xlabel('Degree','FontSize',20);
% ylabel('J_{2n} ','FontSize',20);
% 
% 
% for i=1:size(rhocorei,1)
%     for j=1:size(rhocorei,2)
%         
%         for k=1:size(J,3)
%             Jslice(k)=J(i,j,k);
%         end
%         
%         p=polyfit(2:2:nmax,Jslice,1);
%         
%         pe(i,j)=p(1);
%         
%     end   
% end


% 
% figure; hold on;
% set(gca,'FontSize',20);
% pcolor(rcorei/1000,rhocorei,J(:,:,1)); shading interp;
% colorbar('FontSize',20);
% 
% 
% xlabel('Core size [km] ','FontSize',12);
% ylabel('Core density [kg/m^3] ','FontSize',12);


%% Computing C/MR^2
% 
% [a1,c1]=f2axes(router,fouteri);
% [a2,c2]=f2axes(rcorei,fcorei);
% 
% M1=4/3*pi*rhoouteri*router.^3;
% M2=4/3*pi*(rhocorei-rhoouteri).*rcorei.^3;
% 
% Mt=M1+M2;
% 
% Ch1=0.2*(M1).*(a1.^2+c1.^2);
% Ch2=0.2*(M2).*(a2.^2+c2.^2);
% 
% Ch=Ch1+Ch2;
% 
% lambdah=Ch./(Mt.*router.^2);
% 
% levels=-0.5:0.02:0.6666;
% [C,h]=contour(rcorei/1000,rhocorei,lambdah,levels,'Color','g','LineWidth',2);
% 
% clabel(C,h,'manual','FontSize',fntsize_sm,'Color','k','EdgeColor','k','BackGroundColor','w');
% 
% print(gcf, '-dpsc', 'Ceres2Layer.eps');
% 
% 
% toc







