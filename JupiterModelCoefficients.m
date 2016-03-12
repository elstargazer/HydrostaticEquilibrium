function C = JupiterModelCoefficients(rho_core,r_core,GM_obs,n)

G = 6.67e-11;
R1 = 69911000;
T = 9.925;
K = 1/446;

MaxDegree = 6;

V = 4/3*pi*R1^3;
V_core = 4/3*pi*r_core^3;
rho_mean = GM_obs/G/V;

Rref = R1;
tol = 1;
eps = tol + 1;
%% Calculate equilibrium
M_obs = GM_obs/G;
M_core = 4/3*pi*(r_core.^3)*rho_core;
% F_core = M_core./(M_obs);

M_env = M_obs - M_core;
V_env = V     - V_core;

rho_env_mean = M_env/V_env;

nlayer = 25;
a   = linspace(R1,r_core,nlayer);
rho = [linspace(rho_env_mean,rho_env_mean,nlayer-1)...
    rho_core];

a_midpoint = (a + [a(2:end) 0])/2;

% figure; hold on

while eps > tol
rhod = [rho(1) diff(rho)];
W=2*pi/(3600*T);

Md = 4/3*pi.*rhod.*a.*a.*a;
M  = sum(Md);

mu = a./a(1);
sigma = rhod./rho(1);
Lambda2 = W.*W/(pi*G*rho(1));

A=zeros(nlayer);
D = zeros(nlayer,1);

for i=1:nlayer
    
    D(i) = 20*sum(sigma(1:i-1)) + 8*sigma(i) + ...
        20*sum((mu(i+1:end)/mu(i)).^3.*sigma(i+1:end));
    for j=1:nlayer
        
        if (i==j)
            A(i,j) = 1; 
        elseif (j<i)
            A(i,j) = -12*sigma(j)/D(i);
        elseif (j>i)
            A(i,j) = -12*((mu(j)/mu(i)).^5)*sigma(j)/D(i);
        end    
    end
end

B = 15*Lambda2./D;

e2 = (linsolve(A,B))';
% e = sqrt(e2);
c = a.*sqrt(1-e2);

% figure; hold on;
% plot(a,e);

% fp=1-sqrt(1-e2(1));

% Md = 4/3*pi.*rhod.*a.*a.*c;
% M  = sum(Md);
% GM = M*G;
% g_surf = GM/a(1)/a(1);

%% Compute pressure

% compute pressure at the center and at the boundary of layers

P          = zeros(1,nlayer);
U_top      = zeros(1,nlayer);
U_midpoint = zeros(1,nlayer);

for i = 1:nlayer  
   U_top(i)      = EllPotMultilayer(a,c,a(i),0,0,rhod,W);
   U_midpoint(i) = EllPotMultilayer(a,c,a_midpoint(i),0,0,rhod,W);
end


% x = linspace(-1.5*R1,1.5*R1,20);
% z = x;
% 
% [xi,zi] = meshgrid(x,z);
% 
% Ui = xi*0;
% 
% for l=1:numel(Ui)    
%     Ui(l) = EllPotMultilayer(a,c,xi(l),0,zi(l),rhod,W);
% end
% 
% 
% figure; hold on;
% surf(xi,zi,Ui); shading interp;
% plot3(a,a*0,U_top,'-ok');
% 
% contour(xi/R1,zi/R1,Ui,10);

for i = 1:nlayer   
    P(i) = abs(sum(rhod(1:i).*(U_midpoint(i)-U_top(1:i))));
end

% old way
% P(1) = 0;
% for i = 2:nlayer   
%     P(i) = abs(sum(rhod(1:i-1).*(U_top(i)-U_top(1:i-1))));   
% end

% central pressure
% Uc = EllPotMultilayer(a,c,0,0,0,rhod,W);
% Pc = sum(rhod.*(Uc-U));

%% Force to be a barotrope

rho_pt = EOS_barotrope(P,K,n,a,r_core,rho_core);
rhod_pt = [rho_pt(1) diff(rho_pt)];

% plot(a/1000,rho_pt,'r');

Md_pt = 4/3*pi.*rhod_pt.*a.*a.*c;
M_pt  = sum(Md_pt);
M_env_pt = sum(Md_pt)-M_core;

K = K*(M_env/M_env_pt);

rho_pt = EOS_barotrope(P,K,n,a,r_core,rho_core);
rhod_pt = [rho_pt(1) diff(rho_pt)];

Md_pt = 4/3*pi.*rhod_pt.*a.*a.*c;
M_pt  = sum(Md_pt);

% M_pt/M_obs

eps = max(rho_pt-rho);

rho = rho_pt;
% rho(1) = 1;

% % figure; hold on;
% plot(a,P,'b');

% stairs([a/1000 0],[rho_pt rho_core]);

% stairs([a/1000 0],[rho_pt rho_core]);

% plot(a/1000,rho_pt)

end

%% Plot

% figure; hold on;
% plot([a],[P]/1e9,'o');
% ylabel('Pressure [GPa]');
% 
% figure; hold on;
% plot([a],[rho_pt],'o');
% ylabel('Density [kg/m3]');


%% Compute gravity coefficients

w  = Md_pt./M_obs;

% sum(Md_pt)/M_obs

l = (2:2:MaxDegree)';
l_rep = repmat((2:2:MaxDegree)',[1 size(Md_pt,2)]);
c_rep = repmat(c,[numel(l) 1]);
a_rep = repmat(a,[numel(l) 1]);
w_rep = repmat(w,[numel(l) 1]);

l2 = l_rep/2;

C_rep = 3./(Rref.^(l_rep)) .* ((a_rep.^2-c_rep.^2).^l2) .* ...
    ((-1).^l2) ./ ( (2*l2+1).*(2*l2+3) );   

% C_rep=C_rep.*NormCoef(l_rep,0);
C = sum(w_rep.*C_rep,2);













