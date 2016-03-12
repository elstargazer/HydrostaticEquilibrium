function rho = EOS_barotrope(P,K,n,a,r_core,rho_core)

rho = zeros(size(a));

for i = 1:numel(a)
    
    if a(i) <= r_core + (1e-6)
        rho(i) = rho_core;     
    else
        rho(i) = K*P(i).^(n/(1+n));
    end 
end
