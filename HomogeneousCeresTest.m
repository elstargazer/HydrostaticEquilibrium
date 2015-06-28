ccc;
tic;
%% Input parameters
G=6.67384e-11;

GM=62.68e9; % PCK version 0.2
T=9.073859324514187; % DLR

ell = [484.17 481.41 447.80];
[fap_obs, fq_obs,fabp_obs] = ell2f(ell);
fbp_obs = (ell(2)-ell(3))/ell(2);


fntsize = 12;
fntsize_sm = 8;
im_size=[0 0 13 9];


fig_folder='~/Dawn/Papers/CeresPaper1/';
%% Preliminary computations

M=GM/G;
R=prod(ell).^(1/3);
V=1e9*4/3*pi*R.^3;

% sigma_fouter_obs=sqrt(...
%     (c.*c.*sigmaa.*sigmaa+a.*a.*sigmac.*sigmac)./(c.^4));

router=(V./(4/3*pi)).^(1/3);
rhomean=M./V;
g=GM./(router.^2);

%% Homogeneous Ceres

Th=8:0.1:10;
f0=0.01;

fh_tc=HydrostaticStateExact(router,T,rhomean,f0);

for i=1:numel(Th)
    if (i>1)
        f0=fh(i-1);
    else
        f0=0.001;
    end
    [fh(i),fval]=HydrostaticStateExact(router,Th(i),rhomean,f0);   
end

figure;
set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
hold on;box on;grid on;

phomo=plot(Th,fh,'-k','LineWidth',3);

ylabel('Outer flattening []','FontSize',fntsize);


%% Pre-Dawn shapes

fp = [546, 644, 601, 542, 915, 897, 380, 669, 736]/10000;
fp_err = [103, 62, 550, 454, 184, 174, 158, 45, 62]/10000;

ccj=jet(numel(fp));

for i=1:numel(fp);
    
    plot(8+i/10,fp(i),'o','Color',ccj(i,:));
    p_(i)=plot(NaN,'-','Color',ccj(i,:));
    errorbar(8+i/10,fp(i),fp_err(i),'Color',ccj(i,:));
    
end

plot([T T],get(gca,'YLim'),'-r','LineWidth',2);

% plot(get(gca,'XLim'),[fabp_obs fabp_obs],'--k'); 
% plot(get(gca,'XLim'),[fap_obs fap_obs],'--r'); 
% plot(get(gca,'XLim'),[fbp_obs fbp_obs],'--b'); 
% plot(get(gca,'XLim'), ones(1,2)*(ell(2)-ell(3))/ell(1),'--g'); 
% plot(get(gca,'XLim'), ones(1,2)*(ell(1)-ell(3))/ell(2),'--y'); 


x_fill = [get(gca,'XLim') fliplr(get(gca,'XLim'))];      % repeat x values
y_fill = fliplr([ones(1,2)*(ell(2)-ell(3))/ell(1),...
     ones(1,2)*(ell(1)-ell(3))/ell(2)]);   % vector of upper & lower boundaries
p_fill=fill(x_fill,y_fill,'g','FaceAlpha',0.3,'EdgeColor','none');  



legend([phomo p_ p_fill],...
    {'$f_{p,homo}$',...
    'Millis et al. 1987 (sol 1)',...
    'Millis et al. 1987 (sol 2)',...
    'Saint-Pe et al. 1993',...
    'Mitchell et al. 1996',...
    'Drummond et al. 1998 (sol 1)',...
    'Drummond et al. 1998 (sol 2)',...
    'Parker et al. 2002',...
    'Thomas et al. 2005',...
    'Carry et al. 2008',...
    'Dawn'},...
    'FontSize',fntsize_sm,...
    'interpreter','latex','boxoff')

set(gca,'xaxisLocation','top')
% 
% axes('xlim', get(gca,'XLim'), 'ylim', get(gca,'YLim'),...
%     'color', 'none', 'XAxisLocation', 'bottom','ytick',[],...
%     'FontSize',fntsize);
xlabel('Rotation period [h]','FontSize',fntsize);

PrintWhite([fig_folder 'Fig_CeresHomo.jpg']);







