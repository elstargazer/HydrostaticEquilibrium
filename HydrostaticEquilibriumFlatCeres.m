% close all
ccc;
tic;
%% Input parameters
G    = 6.67384e-11;
GM   = 62.68e9; % PCK version 0.2
T    = 9.073859324514187; % DLR
%  T = 8.46;
Rref = 470000;

Npts = 30;

ell     = [ 483.2541034821433 480.9849831138582  445.9487550311291]*1000;
% ell_eq = [481.85802     481.14698     447.03731]*1000;
% ell     = ell_spc*1000;

ell_err = [0.0 0.0 0.0];
J2obs   = 1.185090505130000e-2;
J4obs   = - 5.7288e-04;

J4hyd = -5.1887e-04


offset  = 0.9694*1000;
% OpNav5 JPL shape

[fpa1_obs, fq1_obs, fp1_obs] = ell2f(ell);
% [fpa1_eq, fq1_eq, fp1_eq] = ell2f(ell_eq);

fpb1_obs = (ell(2)-ell(3))/ell(2);

r2_min = 10000;
r2_max = 470000;
rho2_max = 4000;

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

r2=linspace(r2_min,r2_max,Npts);

r2 = r2(1:end-1);
r2_add = linspace(r2(end),r1,20);
r2 = [r2 r2_add(2:end)];

rho2_add = linspace(rhomean,rhomean+100,20);
rho2=linspace(rhomean+100,rho2_max,Npts);
rho2 = [rho2_add(1:end-1) rho2];

[rho2i,r2i]=meshgrid(rho2,r2);

rho1i=-(3*M-4*pi*(r2i.^3).*rho2i)./(4*pi*(r2i.^3)-4*pi*(r1^3));
rho1i(rho1i<0)=NaN;

M2 = 4/3*pi.*(r2i.^3).*(rho2i-rho1i);
M1 = 4/3*pi.*(r1.^3).*(rho1i);

%% Compute hydrostatic flattening factors
[f2i,f1i]=HydrostaticStateExact2lGrid(r1,r2i,T,rho1i,rho2i);

f1i_n = f1i;
f1i_n((f1i<0) | isnan(rho1i)) = NaN;

a1i = fr2abc(r1,f1i,0);

%% Compute gravity coeffs

[J2hi,J4hi]=RadFlat2J(r1,r2i,f1i,f2i,rho1i,rho2i,Rref);

%% Plotting settings

fig_2l=figure;
set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
hold on;box on;grid on;

xlim([r2_min r2_max]/1000);
ylim([rhomean rho2_max]);

xlabel('Core size [km]','FontSize',fntsize);
ylabel('Core density [kg/m^3]','FontSize',fntsize);

%% Contour mantle density

pcolor(r2i/1000,rho2i,rho1i); shading interp;
cbar=colorbar('FontSize',fntsize);
ylabel(cbar,'Outer density [kg/m^3] ','FontSize',fntsize);
caxis([920 rhomean]);
plot(get(gca,'xlim'),[rhomean rhomean],'--k','LineWidth',3);

%% Contour flattening

% 
% levels=0.045:0.005:0.6;
% [C,h]=contour(r2i/1000,rho2i,f1i_n,levels,'Color',...
%     [0.7 0.0 0.45],'LineWidth',2);
% 
% clabel(C,h,'manual','FontSize',fntsize_sm,'Color','k','EdgeColor','k','BackGroundColor','w');
% 
% error_levels=[fp1_obs-sigma_fp1_obs fp1_obs+sigma_fp1_obs];
% 
% [Chyd,h]=contour(r2i/1000,rho2i,f1i_n,[fp1_obs fp1_obs],...
%     'Color','r','LineWidth',5,'LineStyle','-');
% 
% [Chyd_a,~]=contour(r2i/1000,rho2i,f1i_n,[fpa1_obs fpa1_obs],...
%     'Color','r','LineWidth',5,'LineStyle','-');
% 
% [Chyd_b,~]=contour(r2i/1000,rho2i,f1i_n,[fpb1_obs fpb1_obs],...
%     'Color','r','LineWidth',5,'LineStyle','-');
% 
% [Chyd_err,h]=contour(r2i/1000,rho2i,f1i_n,error_levels,...
%     'Color','r','LineWidth',5,'LineStyle','-');
% 
% level_set1=[920 920];
% level_set2=[1000 1250 1500 1750 2000];
% 
% [C,h]=contour(r2i/1000,rho2i,rho1i,level_set1,'--',...
%     'Color','k','LineWidth',1);
% 
% clabel(C,h,'manual','FontSize',fntsize_sm,...
%     'Color','k','EdgeColor','k','BackGroundColor','w');
% 
% % [C,h]=contour(rcorei/1000,rhocorei,rhoouteri{ModelIn},level_set2,...
% %     'Color','c','LineWidth',4);
% % clabel(C,h,'manual','FontSize',fntsize_sm,'Color','k','EdgeColor','k','BackGroundColor','w');
% 
% cbar=colorbar('FontSize',fntsize);
% ylabel(cbar,'Outer density [kg/m^3] ','FontSize',fntsize);
% 
% PrintWhite(fig_2l,[fig_folder 'Fig_2layer.jpg']);


% legend([pl_shell pl_offset],...
%     {'Shell thickness','Core offset'},...
%     'FontSize',fntsize_sm,'Color','w');

%% Non hydrostatic core

% gi = ginput(1);
%
% rcore_test = gi(1);???
% rho_core_test = gi(2);
% rhoouter_test=-(3*M-4*pi*(rcore_test.^3).*rho_core_test)./...
%     (4*pi*(rcore_test.^3)-4*pi*(router(2)^3));

% [fp_noneq_test,fq_noneq_test,fval,exitflag]=...
%     HydrostaticStateExact2lGridCoreNonEq(...
%     router(2),rcore_test*1000,T,rhoouter_test,...
%     rho_core_test,fpa_obs,fq_obs);

% tic
% [fp_noneq,fq_noneq,fval,exitflag]=...
%     HydrostaticStateExact2lGridCoreNonEq(...
%     r1,r2i,T,rho1i,rho2i,fpa1_obs,fq1_obs);
% toc

% convergence conditions
% conv_cond = log10(fval) < 0;
%
% fq_noneq(~conv_cond)=NaN;
% fp_noneq(~conv_cond)=NaN;
%
% figure;
% set(gcf, 'Units','centimeters', 'Position',im_size)
% set(gcf, 'PaperPositionMode','auto')
% set(gca, 'FontSize',fntsize);
% hold on;box on;grid on;
%
% xlim([10 450]);
% ylim([rhomean 5000]);
% caxis([920 rhomean]);
%
% f_ratio = fp_noneq./fq_noneq;
% f_levels=0:0.1:1;
%
% pcolor(r2i/1000,rho2i,rho1i); shading interp;
%
% contour(r2i/1000,rho2i,fp_noneq,f_levels,...
%     '-','Color','r','ShowText','on');
% h1_ = plot(NaN, '-r');
%
% contour(r2i/1000,rho2i,fq_noneq,f_levels,...
%     '-','Color','b','ShowText','on');
% h2_ = plot(NaN, '-b');
%
% contour(r2i/1000,rho2i,f_ratio,[1 1],'-k',...
%     'LineWidth',3,'HandleVisibility','off');
% h3_ = plot(NaN, '-k','LineWidth',3);
%
% legend([h1_ h2_,h3_],...
%     {'$\Delta f_{p,core}$','$\Delta f_{q,core}$','$f_{q,core}=f_{p,core}$'},...
%     'FontSize',fntsize_sm,'Interpreter','latex');
%
% cbar=colorbar('FontSize',fntsize);
% ylabel(cbar,'Outer density [kg/m^3] ','FontSize',fntsize);
%
% xlabel('Core size [km]','FontSize',fntsize);
% ylabel('Core density [kg/m^3]','FontSize',fntsize);
%
% PrintWhite([fig_folder 'Fig_NonhydroCore.jpg']);

%% Computing J2 and C22 for non-hydrostatic core case
%
% [J2_1,C22_1]=rf2j2c22(r1,fq1_obs,fpa1_obs,Rref);
% [J2_2,C22_2]=rf2j2c22(r2i,fq_noneq,fp_noneq,Rref);
%
% J2 = (J2_2.*M2+J2_1.*M1)./M;
% C22 = (C22_2.*M2+C22_1.*M1)./M;
%
% figure;
% set(gcf, 'Units','centimeters', 'Position',im_size)
% set(gcf, 'PaperPositionMode','auto')
% set(gca, 'FontSize',fntsize);
% hold on; box on;grid on;
%
% xlim([10 450]);
% ylim([rhomean 5000]);
% caxis([920 rhomean]);
%
% pcolor(r2i/1000,rho2i,rho1i); shading interp;
%
% levels = 0:0.00001:max(J2(:));
% contour(r2i/1000,rho2i,J2,levels,'-','Color','r','ShowText','on');
% h1_ = plot(NaN, '-r');
%
% levels = 0:0.00001:max(C22(:));
% contour(r2i/1000,rho2i,C22,levels,'-','Color','b','ShowText','on');
% h2_ = plot(NaN, '-b');
%
% legend([h1_ h2_],...
%     {'$J_{2}$','$C_{22}$'},'FontSize',fntsize_sm,'Interpreter','latex');
%
% PrintWhite([fig_folder 'Fig_GravityPrediction.jpg']);

%% Computing hydrostatic J2 and J4

fig_2l_j2 = figure;
set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
hold on;box on;grid on;

xlim([r2_min r2_max]/1000);
ylim([rhomean rho2_max]);

xlabel('Core size [km]','FontSize',fntsize);
ylabel('Core density [kg/m^3]','FontSize',fntsize);

% plot mantle density colors
pcolor(r2i/1000,rho2i,rho1i); shading interp;

cbar=colorbar('FontSize',fntsize);
ylabel(cbar,'Outer density [kg/m^3] ','FontSize',fntsize);
caxis([920 rhomean]);

plot(get(gca,'xlim'),[rhomean rhomean],'--k','LineWidth',3);

levels=(0:0.001:1);
[C,h]=contour(r2i/1000,rho2i,J2hi,levels,'Color','b','LineStyle','-','LineWidth',1);
clabel(C,h,'manual','FontSize',12,'Color','k','EdgeColor','k','BackGroundColor','w');

% gravity solution
% plot solution from J2
[CJhyd,h]=contour(r2i/1000,rho2i,J2hi,[J2obs J2obs],...
    'Color','b','LineStyle','-','LineWidth',2);

% plot bounds for J2
[CJhyd_plus,~]=contour(r2i/1000,rho2i,J2hi,[J2obs J2obs]*(1+0.0314),...
    'Color','b','LineStyle','--','LineWidth',1);

[CJhyd_minus,~]=contour(r2i/1000,rho2i,J2hi,[J2obs J2obs]*(1-0.0314),...
    'Color','b','LineStyle','--','LineWidth',1);

% shape solution
% plot solution from fp
[Cfphyd,h]=contour(r2i/1000,rho2i,f1i_n,[fp1_obs fp1_obs],...
    'Color','r','LineStyle','-','LineWidth',2);

% plot bounds for fp
[Cfphyd_plus,~]=contour(r2i/1000,rho2i,f1i_n,[fpa1_obs fpa1_obs],...
    'Color','r','LineStyle','--','LineWidth',1);

[Cfphyd_minus,~]=contour(r2i/1000,rho2i,f1i_n,[fpb1_obs fpb1_obs],...
    'Color','r','LineStyle','--','LineWidth',1);

% find gravity solution
ncontours = CJhyd(2,1);
r2_Jh   = CJhyd(1,2:ncontours+1)*1000;
rho2_Jh = CJhyd(2,2:ncontours+1);

ncontours = CJhyd_plus(2,1);
r2_Jh_plus   = CJhyd_plus(1,2:ncontours+1)*1000;
rho2_Jh_plus = CJhyd_plus(2,2:ncontours+1);

ncontours = CJhyd_minus(2,1);
r2_Jh_minus   = CJhyd_minus(1,2:ncontours+1)*1000;
rho2_Jh_minus = CJhyd_minus(2,2:ncontours+1);

rho1_Jh=griddata(r2i,rho2i,rho1i,r2_Jh,rho2_Jh,'linear');
M2_Jh=griddata(r2i,rho2i,M2,r2_Jh,rho2_Jh,'linear');

rho1_Jh_plus=griddata(r2i,rho2i,rho1i,r2_Jh_plus,rho2_Jh_plus,'linear');
M2_Jh_plus=griddata(r2i,rho2i,M2,r2_Jh_plus,rho2_Jh_plus,'linear');

rho1_Jh_minus=griddata(r2i,rho2i,rho1i,r2_Jh_minus,rho2_Jh_minus,'linear');
M2_Jh_minus=griddata(r2i,rho2i,M2,r2_Jh_minus,rho2_Jh_minus,'linear');

% find shape solution

ncontours = Cfphyd(2,1);
r2_fph   = Cfphyd(1,2:ncontours+1)*1000;
rho2_fph = Cfphyd(2,2:ncontours+1);

ncontours = Cfphyd_plus(2,1);
r2_fph_plus   = Cfphyd_plus(1,2:ncontours+1)*1000;
rho2_fph_plus = Cfphyd_plus(2,2:ncontours+1);

ncontours = Cfphyd_minus(2,1);
r2_fph_minus   = Cfphyd_minus(1,2:ncontours+1)*1000;
rho2_fph_minus = Cfphyd_minus(2,2:ncontours+1);

rho1_fph=griddata(r2i,rho2i,rho1i,r2_fph,rho2_fph,'linear');
M2_fph=griddata(r2i,rho2i,M2,r2_fph,rho2_fph,'linear');

rho1_fph_plus=griddata(r2i,rho2i,rho1i,r2_fph_plus,rho2_fph_plus,'linear');
M2_fph_plus=griddata(r2i,rho2i,M2,r2_fph_plus,rho2_fph_plus,'linear');

rho1_fph_minus=griddata(r2i,rho2i,rho1i,r2_fph_minus,rho2_fph_minus,'cubic');
M2_fph_minus=griddata(r2i,rho2i,M2,r2_fph_minus,rho2_fph_minus,'cubic');

%% Computing C/MR^2

[a1,c1]=f2axes(r1,f1i);
[a2,c2]=f2axes(r2i,f2i);

Ch1=0.2*(M1).*(2*a1.^2);
Ch2=0.2*(M2).*(2*a2.^2);

Ch=Ch1+Ch2;

lambdah=Ch./(M.*a1.^2);

fig_2l_solns = figure('Color','w');
set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize,'Color','w',...
    'YColor','k','XColor','k');
hold on;box on;grid on;

xlim([r2_min r2_max]/1000);
ylim([rhomean rho2_max]);

xlabel('Core size [km]','FontSize',fntsize);
ylabel('Core density [kg/m^3]','FontSize',fntsize);

% plot mantle density colors
pcolor(r2i/1000,rho2i,rho1i); shading interp;
cbar=colorbar('FontSize',fntsize,'Color','k');
ylabel(cbar,'Outer density [kg/m^3] ','FontSize',fntsize);
caxis([920 rhomean]);
plot(get(gca,'xlim'),[rhomean rhomean],'--k','LineWidth',3);

% lambda levels
levels=-0.5:0.01:0.6666;
[C_I,h]=contour(r2i/1000,rho2i,lambdah,levels,...
    'Color','g','LineWidth',2);

% shape constraint
contour(r2i/1000,rho2i,f1i_n,[fp1_obs fp1_obs],...
    'Color','r','LineWidth',3,'LineStyle','-');
h_f = plot(NaN,'-r','LineWidth',3);

% plot confidence limits
contour(r2i/1000,rho2i,f1i_n,[fpa1_obs fpa1_obs],...
    'Color','r','LineWidth',1,'LineStyle','--');
h_f_err = plot(NaN,'--r');
contour(r2i/1000,rho2i,f1i_n,[fpb1_obs fpb1_obs],...
    'Color','r','LineWidth',1,'LineStyle','--');
h_f_err = plot(NaN,'--r');

% J2 constraint
contour(r2i/1000,rho2i,J2hi,[J2obs J2obs],...
    'Color','b','LineStyle','-','LineWidth',3);
h_j = plot(NaN,'-b','LineWidth',3);

% plot confidence limits
contour(r2i/1000,rho2i,J2hi,[J2obs J2obs]*(1+0.0314),...
    'Color','b','LineStyle','--','LineWidth',1);
h_j_err = plot(NaN,'-b');

contour(r2i/1000,rho2i,J2hi,[J2obs J2obs]*(1-0.0314),...
    'Color','b','LineStyle','--','LineWidth',1);
h_j_err = plot(NaN,'-b');

%J4 constraint
contour(r2i/1000,rho2i,J4hi,[-5.3334e-04 -5.3334e-04],...
    'Color','y','LineStyle','-','LineWidth',3);
h_j4 = plot(NaN,'-y','LineWidth',2);

% legend([h_f h_j h_js h_Z_sol],{'Shape solution','Gravity solution','Solution for T = 8.46 [h]','Admittance solution'},'FontSize',fntsize,...
%     'Color','w','Color','w','TextColor','k');

clabel(C_I,h,'manual','FontSize',fntsize_sm,...
    'Color','k','EdgeColor','k','BackGroundColor','w');

% PrintWhite(fig_2l_solns,[fig_folder 'Fig_2l_solutions.jpg']);

%% Moment of inertia
lambda_J = griddata(r2i,rho2i,lambdah,r2_Jh,rho2_Jh,'linear');
lambda_f = griddata(r2i,rho2i,lambdah,r2_fph ,rho2_fph,'linear');

fig_lambda=figure;
set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
hold on;grid on; box on;

pl_lambdaJ = plot(r2_Jh/1000,lambda_J,'-','LineWidth',3,'Color','b');
pl_lambdaf = plot(r2_fph/1000,lambda_f,'-','LineWidth',3,'Color','r');

lambda_rd = 0.360;
plot_rd = plot(get(gca,'XLim'),[lambda_rd lambda_rd],'--k');

xlabel('Core size [km]','FontSize',fntsize);
ylabel('\lambda','FontSize',fntsize);

ylim([0.34 0.42]);

legend([pl_lambdaJ pl_lambdaf plot_rd],...
    {'$\lambda_{J_{2}}$','$\lambda_{f_{p}}$','$\lambda_{Radau-Darwin}$'},...
    'FontSize',fntsize_sm,'interpreter','latex');

PrintWhite(fig_lambda,[fig_folder 'Fig_Lambda_k.jpg']);

%% Shell thickness

fig_shell=figure;
set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
hold on;grid on; box on;

h_st_J = plot(rho1_Jh,(r1-r2_Jh)/1000,'-b','LineWidth',3);
plot(rho1_Jh_minus,(r1-r2_Jh_minus)/1000,'--b','LineWidth',1);
plot(rho1_Jh_plus,(r1-r2_Jh_plus)/1000,'--b','LineWidth',1);

h_st_fp = plot(rho1_fph,(r1-r2_fph)/1000,'-r','LineWidth',3);
plot(rho1_fph_minus,(r1-r2_fph_minus)/1000,'--r','LineWidth',1);
plot(rho1_fph_plus,(r1-r2_fph_plus)/1000,'--r','LineWidth',1);

xlim([800 2100]);
ylim([0 150]);

legend([h_st_J h_st_fp],{'Gravity solution','Shape solution'},'FontSize',fntsize_sm); 

xlabel('Shell density [kg/m^3]','FontSize',fntsize);
ylabel('Shell thickness [km]','FontSize',fntsize);

% PrintWhite(fig_shell,[fig_folder 'Fig_IceThickness.jpg']);











