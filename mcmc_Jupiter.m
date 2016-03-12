function param = mcmc_Jupiter(N,s,data,data_std,param0,cost_fun,step_param,GM_obs)

% initialize parameter vectors
param=zeros(N,3);
param(1,:)=param0;

% compute initial cost
cost_c=cost_fun(param(1,:),data,data_std,GM_obs);

% cost_c -- old cost
% cost_p -- new (proposed) cost

for i=1:N
    
    param_p = step_param(param(i,:),s); 
    cost_p = cost_fun_Jupiter(param_p,data,data_std,GM_obs);
    metric=exp(- (cost_p) + (cost_c));
    
    if (rand < metric) && (~isnan(cost_p)) && (~isnan(metric))
        param(i+1,:) = param_p;
        cost_c       = cost_p;
    else
        param(i+1,:) = param(i,:);
    end
end

