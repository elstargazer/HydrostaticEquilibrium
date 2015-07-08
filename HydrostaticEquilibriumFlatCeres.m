% close all
ccc;
tic;
%% Input parameters
G=6.67384e-11;
GM=62.68e9; % PCK version 0.2
T=9.073859324514187; % DLR
Rref=500000;

ell = [484.17 481.41 447.80]*1000;
ell_err = [0.0 0.0 0.0];
offset = 890;
% OpNav5 JPL shape

[fpa1_obs, fq1_obs, fp1_obs] = ell2f(ell);

%% figure settings 

fntsize = 12;
fntsize_sm = 10;
im_size=[0 0 13 9];

fig_folder='~/Dawn/Papers/CeresPaper1/';

%% Preliminary conmutations

M=GM/G;
r1=prod(ell).^(1/3);
% sigma_f1_obs=sqrt(...
%     (c.*c.*sigmaa.*sigmaa+a.*a.*sigmac.*sigmac)./(c.^4));

sigma_fp1_obs = 0;

% 1-outer
% 2-core

V=4/3*pi*r1.^3;

f10=0.1;
f20=0.1;

rhomean=M./V;
g=GM./(r1.^2);

%% Grid of core radii and densities
% rcore=2500:500:470000;
% rhocoreg=2000:100:4000;


r2=linspace(10000,470000,100);
rho2=linspace(rhomean,5000,100);

[rho2i,r2i]=meshgrid(rho2,r2);

rho1i=-(3*M-4*pi*(r2i.^3).*rho2i)./(4*pi*(r2i.^3)-4*pi*(r1^3));
rho1i(rho1i<0)=NaN;


%% Compute hydrostatic flattening factors

[f2i,f1i]=HydrostaticStateExact2lGrid(r1,r2i,T,rho1i,rho2i); 

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

pcolor(r2i/1000,rho2i,rho1i); shading interp;
cbar=colorbar('FontSize',fntsize);
ylabel(cbar,'Outer density [kg/m^3] ','FontSize',fntsize);

caxis([920 rhomean]);
plot(get(gca,'xlim'),[rhomean rhomean],'--k','LineWidth',3);


%% Contour flattening
f1i_n = f1i;
f1i_n((f1i<0) | isnan(rho1i)) = NaN;

levels=0.045:0.005:0.6;
[C,h]=contour(r2i/1000,rho2i,f1i_n,levels,'Color',...
    [0.7 0.0 0.45],'LineWidth',2);

clabel(C,h,'manual','FontSize',fntsize_sm,'Color','k','EdgeColor','k','BackGroundColor','w');

error_levels=[fp1_obs-sigma_fp1_obs fp1_obs+sigma_fp1_obs];

[Chyd,h]=contour(r2i/1000,rho2i,f1i_n,[fp1_obs fp1_obs],...
    'Color','r','LineWidth',5,'LineStyle','-');

[Chyd_err,h]=contour(r2i/1000,rho2i,f1i_n,error_levels,...
    'Color','r','LineWidth',5,'LineStyle','-');

level_set1=[920 920];
level_set2=[1000 1250 1500 1750 2000];

[C,h]=contour(r2i/1000,rho2i,rho1i,level_set1,'--',...
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

M2 = 4/3*pi.*(r2i.^3).*(rho2i-rho1i);
M1 = 4/3*pi.*(r1.^3).*(rho1i);

figure;
set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
hold on;grid on;

ax1 = gca; % current axes
ax1.XColor = 'k';
ax1.YColor = 'k';

r2_h=Chyd(1,:)*1000;
rho2_h=Chyd(2,:);

rho1_h=griddata(r2i,rho2i,rho1i,r2_h,rho2_h,'linear');
M2_h=griddata(r2i,rho2i,M2,r2_h,rho2_h,'linear');

offset_core = M./M2_h*offset;

pl_shell=plot(ax1,rho1_h,(r1-r2_h)/1000,'-','LineWidth',3,'Color','r')

xlabel('Shell density [kg/m^{3}]','FontSize',fntsize);
ylabel('Shell thickness [km]','FontSize',fntsize);
xlim([800 rhomean])


ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none','FontSize',fntsize); 
hold on;

pl_offset=plot(ax2,rho1_h,offset_core/1000,'-','LineWidth',3,'Color','b');

ylabel('Core offset [km]','FontSize',fntsize);
xlim([800 rhomean])

legend([pl_shell pl_offset],...
    {'Shell Thickness [km]','Core offset [km]'},...
    'FontSize',fntsize_sm);

% legend({'SPG','SPC'},'FontSize',fntsize_sm);
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
    r1,r2i,T,rho1i,rho2i,fpa1_obs,fq1_obs);
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

pcolor(r2i/1000,rho2i,rho1i); shading interp;

contour(r2i/1000,rho2i,fp_noneq,f_levels,...
    '-','Color','r','ShowText','on');
h1_ = plot(NaN, '-r');

contour(r2i/1000,rho2i,fq_noneq,f_levels,...
    '-','Color','b','ShowText','on');
h2_ = plot(NaN, '-b');

contour(r2i/1000,rho2i,f_ratio,[1 1],'-k',...
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

[J2_1,C22_1]=rf2j2c22(r1,fq1_obs,fpa1_obs,Rref);
[J2_2,C22_2]=rf2j2c22(r2i,fq_noneq,fp_noneq,Rref);

J2 = (J2_2.*M2+J2_1.*M1)./M;
C22 = (C22_2.*M2+C22_1.*M1)./M;

figure;
set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
hold on; box on;grid on;

xlim([10 450]);
ylim([rhomean 5000]);
caxis([920 rhomean]);

pcolor(r2i/1000,rho2i,rho1i); shading interp;

levels = 0:0.00001:max(J2(:));
contour(r2i/1000,rho2i,J2,levels,'-','Color','r','ShowText','on');
h1_ = plot(NaN, '-r');

levels = 0:0.00001:max(C22(:));
contour(r2i/1000,rho2i,C22,levels,'-','Color','b','ShowText','on');
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







