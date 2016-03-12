function U = EllPotMultilayer(a,c,x,y,z,rhod,W)

nlayer = numel(a);
U = 0;
for i=1:nlayer;
    if (rhod(i) > 1e-6)
    dU = EllPot(a(i),c(i),x,y,z,rhod(i));
    U  = U + dU;   
    end
end

U = U + RotPot(x,y,W);