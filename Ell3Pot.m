function U=Ell3Pot(a,b,c,x,y,z,rho)

w=IsInEllipsoid(a,b,c,x,y,z);

G   = 6.67384e-11;
eps = 1e-12;

if (c>b)
    U=Ell3Pot(a,c,b,x,z,y,rho);
else
    
    if (w<1+eps)
        k=0;
    else
        k=GetK(a,b,c,x,y,z);
    end
    
    % [a a c x y z w]
    
    psi=asin(sqrt((a*a-c*c)/(a*a+k)));
    
    amc2=a*a-c*c;
    amb2=a*a-b*b;
    bmc2=b*b-c*c;
    
    eq2=1-(b/a).^2;
    ep2=1-(c/a).^2;
    
    eratio=eq2/ep2;
    
    [F,E,~] = elliptic12(psi,eratio);
    
    U=pi*G*rho*((2*a*b*c)/sqrt(amc2)*(1-x*x/amb2+y*y/amb2)*F+...
        2*a*b*c/sqrt(amc2)*(x*x/amb2-amc2*y*y/(amb2*bmc2)+z*z/bmc2)*E+...
        2*a*b*c/sqrt((a*a+k)*(b*b+k)*(c*c+k))*((c*c+k)*y*y/bmc2-(b*b+k)*z*z/bmc2));
    
end