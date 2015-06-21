function U=EllPot(a,c,x,y,z,rho)

w=IsInEllipsoid(a,a,c,x,y,z);

G=6.67384e-11;

eps=1e-15;

if (w<1-eps)
    k=0;
elseif ((w>1-eps) & (w<1+eps))
    k=0;
elseif (w>1+eps)
    k=GetK(a,a,c,x,y,z);
else
    k=0;
end

% [a a c x y z w]

psi=asin(sqrt((a*a-c*c)./(a*a+k)));

r=sqrt(x.*x+y.*y+z.*z);

amc2=a*a-c*c;

U=pi*G*rho*((2*a*a*c/sqrt(amc2)).*(1-(r.*r-3*z.*z)./(2*amc2)).*psi+a*a*c/(amc2).*...
    (sqrt(c*c+k)./(a*a+k).*(x.*x+y.*y)-2*z.*z./sqrt(c*c+k)));

