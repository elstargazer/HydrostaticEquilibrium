function [ax, ay, az]=Ell3Acc(a,b,c,x,y,z,rho)

w=IsInEllipsoid(a,b,c,x,y,z);

G=6.67e-11;

eps=1e-12;

if ((w<=1-eps) || (abs(w-1)<eps))
      
    psi=asin(sqrt((a*a-c*c)/(a*a)));
    
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
        2*a*b*c/sqrt((a*a)*(b*b)*(c*c))*((c*c)*2*y/bmc2);
    
    az=2*a*b*c/sqrt(amc2)*(+2*z/bmc2)*E+...
        2*a*b*c/sqrt((a*a)*(b*b)*(c*c))*(-(b*b)*2*z/bmc2);
    
    
else
    
    k=GetK(a,b,c,x,y,z);
    
    psi=asin(sqrt((a*a-c*c)/(a*a + k)));
    
    amc2=a*a-c*c;
    amb2=a*a-b*b;
    bmc2=b*b-c*c;
    
    eq2=1-(b/a).^2;
    ep2=1-(c/a).^2;
    
    eratio=eq2/ep2;
    
    [DkDx,DkDy,DkDz] = GetDK(a,b,c,x,y,z);
    
    [F,E,~] = elliptic12(psi,eratio);
    
    A = (2*a*b*c)/sqrt(amc2)*(1-x*x/amb2+y*y/amb2);
    
    DADx = -4*a.*b.*c.*x./(amb2.*sqrt(amc2));
    DADy =  4*a.*b.*c.*y./(amb2.*sqrt(amc2));
    
%     DFDk = -0.5*amc2./sqrt(amc2.*(a.*a+k).*(c.*c+k).*...
%         (c.*c.*eratio+k-a.*a.*(eratio-1)));
    
    DFDpsi = (1 - eratio .* sin(psi).^2).^(-0.5);
    DpsiDk = -0.5*sqrt(amc2./(c.*c+k))./(a.*a+k);

    DFDk = DFDpsi .* DpsiDk;
    
    DFDx = DFDk .* DkDx;
    DFDy = DFDk .* DkDy;
    DFDz = DFDk .* DkDz;
     
    B = 2*a*b*c/sqrt(amc2)*(x*x/amb2-amc2*y*y/(amb2*bmc2)+z*z/bmc2);
    
    DBDx =  4*a.*b.*c.*x./(amb2.*sqrt(amc2));
    DBDy = -4*a.*b.*c.*sqrt(amc2).*y./(amb2.*sqrt(bmc2));
    DBDz =  4*a.*b.*c.*z./(sqrt(amc2).*bmc2);
   
    DEDpsi = sqrt(1-eratio .* sin(psi).^2);
    
    DEDk = DEDpsi .* DpsiDk;
     
    DEDx = DEDk .* DkDx;
    DEDy = DEDk .* DkDy;
    DEDz = DEDk .* DkDz;
      
    C = 2*a*b*c/sqrt((a*a+k)*(b*b+k)*(c*c+k));    
    D=((c*c+k)*y*y/bmc2-(b*b+k)*z*z/bmc2);
    
    DCDk = -a*b*c*((a*a+k)*(b*b+k)+(a*a+k)*(c*c+k)+(b*b+k)*(c*c+k))./...
        (((a*a+k)*(b*b+k)*(c*c+k)).^(1.5));
    
    DCDx = DCDk * DkDx;
    DCDy = DCDk * DkDy;
    DCDz = DCDk * DkDz;
    
    DDDk =  (y*y-z*z)/bmc2;
    DDDy =  2*(c*c+k)*y/bmc2;
    DDDz = -2*(b*b+k)*z/bmc2;
    
    DDDxtot = DDDk * DkDx;
    DDDytot = DDDk * DkDy + DDDy;
    DDDztot = DDDk * DkDz + DDDz;
    
%     DCDDx1 = DCDx .* D + C .* DDDxtot;
%     DCDDy1 = DCDy .* D + C .* DDDytot;
%     DCDDz1 = DCDz .* D + C .* DDDztot;
    
%     DCDDx = a.*b.*(b+(-1).*c).^(-1).*c.*(b+c).^(-1).*DkDx.*((a.^2+k).*(b.^2+k) ...
%         .*(c.^2+k)).^(-3/2).*(b.^2.*c.^2.*((-1).*c.*y+b.*z).*(c.*y+b.*z)+ ...
%         a.^2.*(b+(-1).*c).*(b+c).*(c.^2.*y.^2+b.^2.*z.^2)+k.*((-2).*c.^4.* ...
%         y.^2+2.*b.^4.*z.^2+b.^2.*c.^2.*((-1).*y.^2+z.^2)+a.^2.*(b+(-1).*c) ...
%         .*(b+c).*(y.^2+z.^2)+k.*((-3).*c.^2.*y.^2+3.*b.^2.*z.^2+k.*((-1).* ...
%         y.^2+z.^2))));
%     
%     DCDDy = a.*b.*(b+(-1).*c).^(-1).*c.*(b+c).^(-1).*((a.^2+k).*(b.^2+k).*( ...
%         c.^2+k)).^(-3/2).*(2.*(a.^2+k).*(b.^2+k).*(c.^2+k).*(2.*(c.^2+k).* ...
%         y+DkDy.*(y+(-1).*z).*(y+z))+DkDy.*(b.^2.*c.^2+a.^2.*(b.^2+c.^2)+ ...
%         2.*(a.^2+b.^2+c.^2).*k+3.*k.^2).*((-1).*c.^2.*y.^2+b.^2.*z.^2+k.*( ...
%         (-1).*y.^2+z.^2))) ;
%     
%     DCDDz = a.*b.*(b+(-1).*c).^(-1).*c.*(b+c).^(-1).*((a.^2+k).*(b.^2+k).*( ...
%         c.^2+k)).^(-3/2).*((-2).*(a.^2+k).*(b.^2+k).*(c.^2+k).*(2.*(b.^2+ ...
%         k).*z+DkDz.*((-1).*y.^2+z.^2))+DkDz.*(b.^2.*c.^2+a.^2.*(b.^2+c.^2) ...
%         +2.*(a.^2+b.^2+c.^2).*k+3.*k.^2).*((-1).*c.^2.*y.^2+b.^2.*z.^2+k.* ...
%         ((-1).*y.^2+z.^2)));

%      ax=DADx .* F + DFDx .* A + DBDx .* E + B .* DEDx + DCDDx;
%      ay=DADy .* F + DFDy .* A + DBDy .* E + B .* DEDy + DCDDy;
%      az=0         + DFDz .* A + DBDz .* E + B .* DEDz + DCDDz;

     ax=DADx .* F + DFDx .* A + DBDx .* E + B .* DEDx + DCDx .* D + C .* DDDxtot;
     ay=DADy .* F + DFDy .* A + DBDy .* E + B .* DEDy + DCDy .* D + C .* DDDytot;
     az=0         + DFDz .* A + DBDz .* E + B .* DEDz + DCDz .* D + C .* DDDztot;   
end
% ax = pi*G*rho*ax;
% ay = pi*G*rho*ay;
% az = pi*G*rho*az;