% close all
ccc;
tic;
%% Input parameters
G=6.67384e-11;

GM=7.210443e9;
T=1.370218*24;

Rref=252100;

fntsize=12;
fntsize_sm=10;
imsize=[13 9];

%% Preliminary conmutations

M=GM/G;
% R=(a.*a.*c).^(1/3);
R=Rref;
V=4/3*pi*R.^3;

% fouter_obs=(a-c)./a;
% sigma_fouter_obs=sqrt(...
%     (c.*c.*sigmaa.*sigmaa+a.*a.*sigmac.*sigmac)./(c.^4));

router=(V./(4/3*pi)).^(1/3);
rhomean=M./V;

%% Grid of core radii and densities
Npts = 30;

rcore=linspace(1000,Rref-1000,Npts);
rhocore=linspace(rhomean,3000,Npts);

[rhocorei,rcorei]=meshgrid(rhocore,rcore);

rhoouteri=-(3*M-4*pi*(rcorei.^3).*rhocorei)./(4*pi*(rcorei.^3)-4*pi*(router^3));
rhoouteri(rhoouteri<0)=NaN;

order = 4; 

[fp1a,fp2a,fq1a,fq2a]=HydrostaticStateAn2lTidalGrid...
    (router,rcorei,T,rhoouteri,rhocorei,order);

tic
[fp1n,fp2n,fq1n,fq2n]=HydrostaticStateExact2lTidalGrid...
    (router,rcorei,T,rhoouteri,rhocorei);
toc

%% Compute gravity coeffs

[J2a,C22a]=RadFlat2J2Tri(router,rcorei,...
    fp1a,fp2a,fq1a,fq2a,rhoouteri,rhocorei,Rref);

[J2n,C22n]=RadFlat2J2Tri(router,rcorei,...
    fp1n,fp2n,fq1n,fq2n,rhoouteri,rhocorei,Rref);

fsa=fp1a./(fp1a-fq1a);
fsn=fp1n./(fp1n-fq1n);

fga=(J2a/NormCoef(2,0))./(C22a/NormCoef(2,2));
fgn=(J2n/NormCoef(2,0))./(C22n/NormCoef(2,2));

%% Plotting 

fig1=figure; hold on;
ax1=gca;
set(ax1,'FontSize',fntsize);
box on;

xlim([0 R]/1000);
ylim([rhomean 3000]);

xlabel('Core size [km]','FontSize',fntsize);
ylabel('Core density [kg/m^3]','FontSize',fntsize);

% plot shell density
pcolor(rcorei/1000,rhocorei,rhoouteri); shading interp;
cbar=colorbar('FontSize',fntsize);
ylabel(cbar,'Outer density [kg/m^3] ','FontSize',fntsize);
caxis([920 rhomean]);

% plot homogeneous line
plot(get(gca,'xlim'),[rhomean rhomean],'--k','LineWidth',3);
set(gcf, 'Units','centimeters', 'Position',[0 0 imsize])

% plot hydrostatic ratio

% levels=(1:0.01:8);
% contour(rcorei/1000,rhocorei,real(fga),levels,...
%     'Color','k','ShowText','on');

levels=(1:0.001:8);
contour(rcorei/1000,rhocorei,real(fgn),levels,...
    'Color','m','ShowText','on');

% pcolor(rcorei/1000,rhocorei,real(fga)); shading interp


