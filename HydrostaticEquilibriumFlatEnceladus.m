% close all
ccc;
tic;
%% Input parameters
G=6.67384e-11;

GM=2.5026e9;
T=0.942*24;

Rref=198200;

fntsize=20;
fntsize_sm=18;
imsize=[54 38];

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
rcore=20000:7500:Rref;
rhocoreg=rhomean:100:3000;

[rhocorei,rcorei]=meshgrid(rhocoreg,rcore);

rhoouteri=-(3*M-4*pi*(rcorei.^3).*rhocorei)./(4*pi*(rcorei.^3)-4*pi*(router^3));
rhoouteri(rhoouteri<0)=NaN;

M=4;

[fp1,fp2,fq1,fq2]=HydrostaticStateAn2lTidalGrid...
    (router,rcorei,T,rhoouteri,rhocorei,M);

tic
[fp1n,fp2n,fq1n,fq2n]=HydrostaticStateExact2lTidalGrid...
    (router,rcorei,T,rhoouteri,rhocorei);
toc

% [fp1,fp2,fq1,fq2]=HydrostaticStateAn2lTidalGrid...
%     (router,100000,T,1000,2000,4);

%% Plotting settings
fig1=figure; hold on;
ax1=gca;
set(ax1,'FontSize',fntsize);
box on;

xlim([10 250]);
ylim([1000 3000]);

xlabel('Core size [km]','FontSize',fntsize);
ylabel('Core density [kg/m^3]','FontSize',fntsize);

%% Contour mantle density
%
pcolor(rcorei/1000,rhocorei,rhoouteri); shading interp;
cbar=colorbar('FontSize',fntsize);
ylabel(cbar,'Outer density [kg/m^3] ','FontSize',fntsize);

caxis([920 rhomean]);
plot(get(gca,'xlim'),[rhomean rhomean],'--k','LineWidth',3);
set(gcf, 'Units','centimeters', 'Position',[0 0 imsize])

%%

fs=fp1./(fp1-fq1);
fsn=fp1n./(fp1n-fq1n);

% levels=0:0.001:10;
% contour(rcorei/1000,rhocorei,real(fp2),levels,...
% 'Color','k','ShowText','on');
%
levels=1:0.01:8;
% contour(rcorei/1000,rhocorei,(fs),levels,...
%     'Color','k','ShowText','on');
%
contour(rcorei/1000,rhocorei,(fsn),levels,...
    'Color','k','ShowText','on','Color','m');

%%

[J2,C22]=RadFlat2J2Tri(router,rcorei,...
    fp1,fp2,fq1,fq2,rhoouteri,rhocorei,Rref);

[J2n,C22n]=RadFlat2J2Tri(router,rcorei,...
    fp1n,fp2n,fq1n,fq2n,rhoouteri,rhocorei,Rref);

fg=J2./C22;
fgn=J2n./C22n;

levels=-(1:0.01:8);
contour(rcorei/1000,rhocorei,real(fg),levels,...
    'Color','k','ShowText','on');

contour(rcorei/1000,rhocorei,real(fgn),levels,...
    'Color','m','ShowText','on');


