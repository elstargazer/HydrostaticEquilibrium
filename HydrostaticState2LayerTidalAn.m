function fh=HydrostaticState2LayerTidalAn(r1,r2,T,rho1,rho2,Order)

if (isnan(rho1) || isnan(rho2) || (rho1<0) || (rho2<0))
    fh=[NaN NaN NaN NaN];
    return
end

tol_in=1e-9;
tol_out=1e-4;

itermax_in=100;
itermax_out=100;

G=6.67e-11;

% 1 - outer
% 2 - inner

Omega=2*pi./(T*3600);
Lambda2=(Omega.^2)./(pi*G*rho1);

sigma2=(rho2-rho1)/rho1;

switch Order
    
    %% Zeroth Order
    
    case 0
        iter_out=0;
        eps_out=tol_out+1;
       
        a2=r2;
        a1=r1;
        
        while (eps_out>tol_out)
            
            mu2=a2/a1;
            
            ep1s=60*Lambda2./(8+20*mu2.^3.*sigma2);
            ep2s=60*Lambda2./(20+8*sigma2);
            
            eq1s=45*Lambda2./(8+20*mu2.^3.*sigma2);
            eq2s=45*Lambda2./(20+8*sigma2);
            
            fp1=1-sqrt(1-ep1s);
            fp2=1-sqrt(1-ep2s);    
            fq1=1-sqrt(1-eq1s);
            fq2=1-sqrt(1-eq2s);
            
            b1=a1.*(1-fq1);
            c1=a1.*(1-fp1); 
            b2=a2.*(1-fq2);
            c2=a2.*(1-fp2);
            
            r1n=(a1.*b1.*c1).^(1/3);
            r2n=(a2.*b2.*c2).^(1/3);
            
            delta_r1=r1-r1n;
            delta_r2=r2-r2n;
            
            a1=a1+delta_r1;
            a2=a2+delta_r2;
            
            eps_out=max(abs([delta_r1 delta_r2]));
            
            iter_out=iter_out+1;
            
            if (iter_out>itermax_out);
                break
            end 
        end
        
        fh(1)=fp1;
        fh(2)=fp2;
        fh(3)=fq1;
        fh(4)=fq2;
        
%% Second  Order
        
    case 2
                
        a2=r2;
        a1=r1;
        
        eps_out=tol_out+1;
        iter_out=0;
        
        while (eps_out>tol_out)
            
            mu2=a2/a1;
            
            Fp=1+0.4*sigma2+2.5*mu2.^3.*sigma2+mu2.^3.*sigma2.^2-...
                (9/10)*mu2.^5.*sigma2;
           
            ep1s=7.5*Lambda2.*(1+0.4*sigma2+0.6*mu2.^5.*sigma2)./Fp;
            ep2s=7.5*Lambda2.*(1+mu2.^3.*sigma2)./Fp;
            
            eq1s=(45/8)*Lambda2.*(1+0.4*sigma2+0.6*mu2.^5.*sigma2)./Fp;
            eq2s=(45/8)*Lambda2.*(1+mu2.^3.*sigma2)./Fp;
            
%             J2=(15/16)*Lambda2*(1+0.4*sigma2+...
%                 (8/5)*mu2.^5.*sigma2+mu2.^8.*sigma2.^2)./...
%                 ((1+mu2.^3.*sigma2).*Fp);
%             
%             C22=(9/32)*Lambda2.*(1+0.4*sigma2+...
%                 (8/5)*mu2.^5.*sigma2+mu2.^8.*sigma2.^2)./...
%                 ((1+mu2.^3.*sigma2).*Fp);
            
            fp1=1-sqrt(1-ep1s);
            fp2=1-sqrt(1-ep2s);  
            fq1=1-sqrt(1-eq1s);
            fq2=1-sqrt(1-eq2s);
            
            b1=a1.*(1-fq1);
            c1=a1.*(1-fp1);  
            b2=a2.*(1-fq2);
            c2=a2.*(1-fp2);
            
            r1n=(a1.*b1.*c1).^(1/3);
            r2n=(a2.*b2.*c2).^(1/3);
            
            delta_r1=r1-r1n;
            delta_r2=r2-r2n;
            
            a1=a1+delta_r1;
            a2=a2+delta_r2;
            
            eps_out=max(abs([delta_r1 delta_r2]));
            
            iter_out=iter_out+1;
            
            if (iter_out>itermax_out);
                break
            end  
        end
        
        fh(1)=fp1;
        fh(2)=fp2;
        fh(3)=fq1;
        fh(4)=fq2;
%% Forth  Order

    case 4
                
        a2=r2;
        a1=r1;
        
        eps_out=tol_out+1;
        iter_out=0;
        
        while (eps_out>tol_out)
            
            fh0=HydrostaticState2LayerTidalAn(r1,r2,T,rho1,rho2,2);
            
            mu2=a2/a1;
%             ep1s=60*Lambda2./(8+20*mu2.^3.*sigma2);
%             ep2s=60*Lambda2./(20+8*sigma2);
%             eq1s=45*Lambda2./(8+20*mu2.^3.*sigma2);
%             eq2s=45*Lambda2./(20+8*sigma2);
            
            fp1=fh0(1);
            fp2=fh0(2);
            fq1=fh0(3);
            fq2=fh0(4);
            
            ep1s=2*fp1-fp1.*fp1;
            ep2s=2*fp2-fp2.*fp2;
            eq1s=2*fq1-fq1.*fq1;
            eq2s=2*fq2-fq2.*fq2;
              
            eps_inn=tol_in+1;
            iter_in=0;
            while (eps_inn>tol_in)
                
                ep1s_old=ep1s;
                ep2s_old=ep2s;
                eq1s_old=eq1s;
                eq2s_old=eq2s;
                
                ep1s=(60*Lambda2+12*mu2.^5.*sigma2.*ep2s_old)./...
                    (8+20*mu2.^3.*sigma2)-...
                    ((22+140.*mu2.^3.*sigma2).*ep1s_old.^2-...
                    (70*mu2.^3+105*mu2.^5).*sigma2.*ep1s_old.*ep2s_old+...
                    (42.*mu2.^5+15.*mu2.^7).*sigma2.*ep2s_old.^2+...
                    16.*ep1s.*eq1s+...
                    (70*mu2.^3-42*mu2.^5).*sigma2.*ep1s_old.*eq2s_old-...
                    (42*mu2.^5-30*mu2.^7).*sigma2.*ep2s_old.*eq2s_old)./...
                    (7*(8+20*mu2.^3.*sigma2));
                                
                ep2s=(60*Lambda2+12*ep1s_old)./(20+8*sigma2)+...
                    (48*ep1s_old.^2-35*ep1s_old.*ep2s_old-...
                    (35+22.*sigma2).*ep2s_old.^2-...
                    12*ep1s_old.*eq1s_old-28*eq1s_old.*ep2s_old-...
                    16*sigma2.*ep2s_old.*eq2s_old)./...
                    (7*(20+8*sigma2));
                
                eq1s=(45*Lambda2+12*mu2.^5.*sigma2.*eq2s_old)./...
                    (8+20*mu2.^3.*sigma2)-...
                    ((8+105*mu2.^3.*sigma2).*eq1s_old.^2+...
                    (42*mu2.^5+15.*mu2.^7).*sigma2.*eq2s_old.^2+...
                    (16.*ep1s_old.*eq1s_old+...
                    (70.*mu2.^3-42*mu2.^5).*sigma2.*eq1s_old.*ep2s_old+...
                    (70*mu2.^3+84*mu2.^5).*sigma2.*eq1s_old.*eq2s_old-...
                    (42*mu2.^5-30*mu2.^7).*sigma2.*ep2s_old.*eq2s_old))./...
                    (7*(8+20*mu2.^3.*sigma2));
                      
                eq2s=(45*Lambda2+12.*eq1s_old)/(20+8*sigma2)+...
                    (48*eq1s_old.^2-8*sigma2.*eq2s_old.^2-...
                    12*ep1s_old.*eq1s_old+...
                    56*eq1s_old.*eq2s_old-28*ep1s_old.*eq2s_old-...
                    16*sigma2.*ep2s_old.*eq2s_old)./...
                    (7*(20+8*sigma2));
                
                eps_inn=max(abs([ep1s_old - ep1s ...
                    ep2s_old - ep2s ...
                    eq1s_old - eq1s ...
                    eq2s_old - eq2s]));
                
                iter_in=iter_in+1;
                
                if (iter_in>itermax_in);
                    break
                end  
            end
            
            fp1=1-sqrt(1-ep1s);
            fp2=1-sqrt(1-ep2s);   
            fq1=1-sqrt(1-eq1s);
            fq2=1-sqrt(1-eq2s);
            
            b1=a1.*(1-fq1);
            c1=a1.*(1-fp1);       
            b2=a2.*(1-fq2);
            c2=a2.*(1-fp2);
            
            r1n=(a1.*b1.*c1).^(1/3);
            r2n=(a2.*b2.*c2).^(1/3);
            
            delta_r1=r1-r1n;
            delta_r2=r2-r2n;
            
            a1=a1+delta_r1;
            a2=a2+delta_r2;
            
            eps_out=max(abs([delta_r1 delta_r2]));
            
            iter_out=iter_out+1;
            
            if (iter_out>itermax_out);
                disp('read max iterations');
                break
            end  
        end
        
        fh(1)=fp1;
        fh(2)=fp2;
        fh(3)=fq1;
        fh(4)=fq2;
        
    otherwise
        error('Unsupported order');      
end

