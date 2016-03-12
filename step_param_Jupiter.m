function param_new=step_param_Jupiter(param,s)

% rho_core=param(1);
% r_core=param(2);
% n = param(3);

param_new = param + s.*randn(size(param));

