function d2=DeltaSquaredRing(f,r,T,rho,Rr,Mr)

W=2*pi/(3600*T);

a=r/((1-f).^(1/3));
c=r/((1-f).^(1/3))-f*r/((1-f).^(1/3));

% fi=(-89:1:89)/180*pi;
% 
% [xa,ya,za]=sph2cart(0,fi,1);
% 
% xa=xa*a;
% ya=ya*a;
% za=za*c;
% 
% Ua=zeros(1,numel(xa));
% 
% for i=1:numel(xa)
%     Ua(i)=EllPotTotRing(a,c,xa(i),ya(i),za(i),rho,W,Rr,Mr);
% end
% 
% d2=max(Ua)-min(Ua);

d2=(EllPotTotRing(a,c,0,0,c,rho,W,Rr,Mr)-EllPotTotRing(a,c,a,0,0,rho,W,Rr,Mr)).^2;