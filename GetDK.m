function [kx,ky,kz]=GetDK(a,b,c,x,y,z)

k=GetK(a,b,c,x,y,z);

bot= ((x./(a.*a+k)).^2)+...
     ((y./(b.*b+k)).^2)+...
     ((z./(c.*c+k)).^2);
           
kx = 2*x./(bot.*(a.*a+k));
ky = 2*y./(bot.*(b.*b+k));
kz = 2*z./(bot.*(c.*c+k));


% kxold = 2*x;
% kyold = 2*y;
% kzold = 2*z;
% 
% tol = 1e-8;
% eps = tol+1;
% 
% 
% while (eps>tol)
%     
%     kxnew = 2*x*k/(a*a+k)+(a*a*x*x/((a*a+k).^2)+...
%                            b*b*y*y/((b*b+k).^2)+...
%                            c*c*z*z/((c*c+k).^2)).*kxold;
%                        
%     kynew = 2*y*k/(b*b+k)+(a*a*x*x/((a*a+k).^2)+...
%                            b*b*y*y/((b*b+k).^2)+...
%                            c*c*z*z/((c*c+k).^2)).*kyold;
%                        
%                        
%     kznew = 2*z*k/(c*c+k)+(a*a*x*x/((a*a+k).^2)+...
%                            b*b*y*y/((b*b+k).^2)+...
%                            c*c*z*z/((c*c+k).^2)).*kzold;
%                        
%                        
%     eps = max(abs([kxold-kxnew kyold-kynew kzold-kznew]));
%     
%     kxold = kxnew;
%     kyold = kynew;
%     kzold = kznew;
%     
% end
% 
% kx = kxnew;
% ky = kynew;
% kz = kznew;
% 
% delta = .01;
% 
% kxp=GetK(a,b,c,x+delta,y,z);
% kxm=GetK(a,b,c,x-delta,y,z);
% 
% kyp=GetK(a,b,c,x,y+delta,z);
% kym=GetK(a,b,c,x,y-delta,z);
% 
% kzp=GetK(a,b,c,x,y,z+delta);
% kzm=GetK(a,b,c,x,y,z-delta);
% 
% kx = (kxp-kxm)/(2*delta);
% ky = (kyp-kym)/(2*delta);
% kz = (kzp-kzm)/(2*delta);




