function [ax, ay, az]=Acc3Pot(a,b,c,x,y,z,rho)

w=IsInEllipsoid(a,b,c,x,y,z);

G=6.67e-11;

eps=1e-12;

if (w<1-eps)
    k=0;
elseif (abs(w-1)<eps)
    k=0;
elseif (w>1+eps)
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

ax=(2*a*b*c)/sqrt(amc2)*(-2*x/amb2)*F+...
    2*a*b*c/sqrt(amc2)*(2*x/amb2)*E;

ay=(2*a*b*c)/sqrt(amc2)*(+2*y/amb2)*F+...
    2*a*b*c/sqrt(amc2)*(-amc2*2*y/(amb2*bmc2))*E+...
    2*a*b*c/sqrt((a*a+k)*(b*b+k)*(c*c+k))*((c*c+k)*2*y/bmc2);

az=2*a*b*c/sqrt(amc2)*(+2*z/bmc2)*E+...
    2*a*b*c/sqrt((a*a+k)*(b*b+k)*(c*c+k))*(-(b*b+k)*2*z/bmc2);

% ax = pi*G*rho*ax;
% ay = pi*G*rho*ay;
% az = pi*G*rho*az;