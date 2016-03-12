function cost=cost_fun_Jupiter(param,C_obs,C_obs_std,GM_obs)

% rho_core=param(1);
% r_core=param(2);
% n = param(3);

C_test = JupiterModelCoefficients(param(1),param(2),GM_obs,param(3));

% cost=((f(2) - f_obs).^2)./(0.004^2) + ...
%      ((J2  - J2_obs).^2)./(25e-5^2);

cost = sum(((C_test - C_obs).^2)./(C_obs_std.^2));


