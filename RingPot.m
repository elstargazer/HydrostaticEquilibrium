function U=RingPot(x,y,z,R,M)

G=6.67e-11;

r=sqrt(x.*x+y.*y);

% rbar=r./R;
% zbar=z./R;
% A=(rbar-1).^2+(z-zbar).^2;
% k=4*rbar./A;
% [K,~] = ellipke((k),eps(1));
% U=-2*G.*M./(pi*R.*sqrt(A)).*K;



q2=R.*R+2*r*R+r.*r+z.*z;

k=(4*R*r./(q2));

[K,~] = ellipke((k),1e-10);

U=2*G*M./(pi*sqrt(q2)).*K;
