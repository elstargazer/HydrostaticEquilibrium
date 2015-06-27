ccc;
tic;
%% Input parameters
G=6.67384e-11;

GM=62.68e9; % PCK version 0.2
T=9.073859324514187; % DLR

ell = [484.17 481.41 447.80];
[fp_obs, fq_obs] = ell2f(ell);

fntsize = 12;
fntsize_sm = 10;
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

plot(Th,fh,'-k','LineWidth',3);
plot([T T],get(gca,'YLim'),'-r','LineWidth',2);

xlabel('Rotation period [h]','FontSize',fntsize);
ylabel('Outer flattening []','FontSize',fntsize);
legend({'f_{p,homo} ','rotation rate'},'FontSize',fntsize_sm)
PrintWhite([fig_folder 'Fig1.eps']);


