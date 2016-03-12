ccc

G = 6.67e-11;
%% general parameters

R1 = 69911000;
R_Earth = 6371000;
M_Earth =  5.9736e24;

N  = 65000;
Nc = 8;

%% Define test Juputer
r_core = 1.5*R_Earth;
rho_core = 5000;
rho_mean = 1326;
GM_obs = 126.687e15;
poly_index = 1;

M_core_obs = 4/3*pi*(r_core.^3).*rho_core;

%% Define gravity field uncertainty

C_obs_std = [1e-6; 1e-6; 1e-6];

%% Compute test case gravity field
tic
C_true = JupiterModelCoefficients(rho_core,r_core,GM_obs,poly_index);
toc

%% Perturb gravity field
C_obs = C_true + 0*randn(size(C_true)).*C_obs_std;

%% Run MCMC to recover test case 

param_start = [rho_core r_core poly_index];
step_sigma  = [400 200000 0.01];

tic
parfor j=1:Nc
    
    param_start_now = param_start + randn(size(step_sigma)).*1.*step_sigma;
    param{j}    = mcmc_Jupiter(N, step_sigma, C_obs, C_obs_std, param_start_now, ...
        @cost_fun_Jupiter,@step_param_Jupiter,GM_obs);  
    j
    
end
toc

% collects all the chains
param_all=[];
for i=1:numel(param)    
    param_all=[param_all; param{i}];    
end


%% Plot heat maps


% box_cond = (param_all(:,1) < min_core_rad) | ...
%            (param_all(:,1) > max_core_rad) | ...
%            (param_all(:,2) < min_core_rho) | ...
%            (param_all(:,2) > max_core_rho);
%        
% param_all(box_cond,:) = [];       

% n=100;
% xi = linspace(min_core_rad,max_core_rad,n);
% yi = linspace(min_core_rho,max_core_rho,n);
% 
% % xi = linspace(min(param_all(:,1)),max(param_all(:,1)),n);
% % yi = linspace(min(param_all(:,2)),max(param_all(:,2)),n);
% 
% xr = interp1(xi,1:numel(xi),param_all(:,1),'nearest')';
% yr = interp1(yi,1:numel(yi),param_all(:,2),'nearest')';
% 
% z = accumarray([xr' yr'], 1, [n n]);
% 
% [xii,yii]=meshgrid(xi,yi);
% 
% figure
% xlim([min_core_rad max_core_rad]/1000);
% ylim([min_core_rho max_core_rho]);
% set(gca,'FontSize',12);
% box on
% hold on
% 
% pcolor(xii/1000,yii,z'); shading interp;


%% 


figure; hold on;
hist(param_all(:,1),30);
xlabel('Core density kg/m^{3}');

figure; hold on;
hist(param_all(:,2)/R_Earth,30);
xlabel('Core radius [R_{Earth}]');

figure; hold on;
hist(param_all(:,3),30);
xlabel('Polytrope index');

rho_core_all = param_all(:,1);
r_core_all   = param_all(:,2);
n_all        = param_all(:,3);
M_core_all = 4/3*pi*(r_core_all.^3).*rho_core_all;


figure; hold on;
hist(M_core_all/M_Earth,30);
xlabel('Core mass [M_{Earth}]');

% Heat map 

min_core_mass = min(M_core_all);
max_core_mass = max(M_core_all);

min_n = min(n_all);
max_n = max(n_all);

n=20;
xi = linspace(min_core_mass,max_core_mass,n);
yi = linspace(min_n,max_n,n);

xr = interp1(xi,1:numel(xi),M_core_all,'nearest')';
yr = interp1(yi,1:numel(yi),n_all,'nearest')';

z = accumarray([xr' yr'], 1, [n n]);

[xii,yii]=meshgrid(xi,yi);

figure
xlim([min_core_mass max_core_mass]/M_Earth);
ylim([min_n max_n]);
set(gca,'FontSize',12);
box on
hold on

pcolor(xii/M_Earth,yii,z'); shading interp;
plot(M_core_obs/M_Earth,poly_index,'ow','MarkerSize',10);

xlabel('Core mass [M_{Earth}]','FontSize',12);
ylabel('Polytropic index','FontSize',12);

%%
% 
% % plot coefficinets
% figure; hold on;
% set(gca,'YScale','log');
% 
% l = 2:2:2*numel(C_obs)
% plot(l,abs(C),'-ok','LineWidth',3,'MarkerSize',5);
% 
% xlabel('Degree');
% ylabel('Coefficient value');


