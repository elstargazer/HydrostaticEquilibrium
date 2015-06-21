function fh=HydrostaticState2LayerAn(r1,r2,T,rho1,rho2,f10,f20,Order)

if (isnan(rho1) || isnan(rho2) || (rho1<0) || (rho2<0))
    fh=[NaN NaN];
    return
end
    

G=6.67e-11;

% 1 - core
% 2 - outer

Omega=2*pi./(T*3600);

V2=4/3*pi*r2.^3;
V1=4/3*pi*r1.^3;

Vsil=V2-V1;
Msil=rho1*Vsil;
Mcore=rho1*V1;

Mtotal=Mcore+Msil;
rhomean=Mtotal/V2;

switch Order
    
    case 1
        
        deltaV=(r1/r2).^3.*(1-rho2/rho1);
        gammaV=2/5+3/5*rho2/rho1;
        Hhydro=(2/5*rhomean/rho1)*(1-3/5*deltaV/gammaV*(r1/r2).^2)/...
            (deltaV+2/5*rho2/rho1-9/25*deltaV*rho2/(gammaV*rho1)*(r1/r2).^2);
        Khydro=(Hhydro-2/5)/(3/5);
        betaV=Omega.*Omega./(pi*G*rhomean);
        fh(1,:)=15/16*Hhydro*betaV*(1-5/16*betaV);
        T2V=Hhydro*5/2*betaV*3/4;
        constA=5/2*(rhomean/rho1)*3/4*betaV;
        constC=3/2*rho2/rho1;
        S2V=(constA+constC*T2V)./(1+constC);
        fh(2,:)=S2V./T2V.*fh(1,:);
        
    case 2
        
        sigma2=(rho1-rho2)/rho2;
        
        a2=r2;
        a1=r1;
        
        eps=1;
        tol=0.001;
        
        Lambda2=(Omega.^2)./(pi*G*rho2);
        
        while (eps>tol)
            
            mu2=a1/a2;
            
            Fp=1+0.4*sigma2+2.5*(mu2.^3).*sigma2+(mu2.^3).*(sigma2.^2)-...
                0.9.*(mu2.^5).*sigma2;
            
            
            e2s=1.875*Lambda2.*(1+0.4.*sigma2+0.6.*mu2.^5.*sigma2)./Fp;
            
            e1s=1.875*Lambda2.*(1+mu2.^3.*sigma2.^2)./Fp;
            
            f1=1-sqrt(1-e1s);
            f2=1-sqrt(1-e2s);
            
            
            b1=a1.*(1-f1);
            b2=a2.*(1-f2);
            
            
            r1n=(a1.*a1.*b1).^(1/3);
            
            r2n=(a2.*a2.*b2).^(1/3);
            
            
            delta_r1=r1-r1n;
            delta_r2=r2-r2n;
            
            a1=a1+delta_r1;
            a2=a2+delta_r2;
            
            eps=max(abs([delta_r1 delta_r2]));
               
        end
        
        fh(2,:)=f2;
        fh(1,:)=f1;
        
    case 6
        
        fh0=HydrostaticState2LayerAn(r1,r2,T,rho1,rho2,f10,f20,2);
        
        f10=fh0(1,:);
        f20=fh0(2,:);
        
        a1=r1./((1-f10).^(1/3));
        a2=r2./((1-f20).^(1/3));
        
        sigma2=(rho1-rho2)/rho2;
        
        a2=r2;
        a1=r1;
        
        Lambda2=(Omega.^2)./(pi*G*rho2);
        
        e1s=2.*f10-f10.^2;
        e2s=2.*f20-f20.^2;
        
        tol1=1e-8;
        tol2=1e-8;
        
        eps2=1;
        i2=1;
        
        while (eps2>tol2)
            
            mu2=a1/a2;
            
            eps1=1;
            i1=1;
            
            while (eps1>tol1)
                
                e1s0=e1s;
                e2s0=e2s;
                
                e2s=(15*Lambda2+12*mu2.^5.*sigma2.*e1s)./(8+20*mu2.^3.*sigma2)-...
                    ((8+105*mu2.^3.*sigma2).*e2s.^2-(70.*mu2.^3+84*mu2.^5).*sigma2.*e2s.*e1s+(42*mu2.^5+15*mu2.^7).*sigma2.*e1s.^2)./...
                    (7*(8+20*mu2.^3.*sigma2))-...
                    1./(14*(8+20*mu2.^3.*sigma2)).*...
                    (175*mu2.^3.*sigma2.*e2s.^3-(105*mu2.^3+210*mu2.^5).*sigma2.*e2s.^2.*e1s-(35*mu2.^3-84*mu2.^5-120*mu2.^7).*sigma2.*e2s.*e1s.^2+...
                    (21*mu2.^5-15*mu2.^7-35*mu2.^9).*sigma2.*e1s.^3);
                
                e1s=(15*Lambda2+12*e2s)./(20+8*sigma2)+(48*e2s.^2-56*e1s.*e2s-8*sigma2.*e1s.^2)./(7*(20+8*sigma2))+...
                    32*(e2s.^3)./(7*(20+8*sigma2));
                
                eps1=max(abs([e1s0-e1s e2s0-e2s]));
                i1=i1+1;
                if (i1>10)
                    break
                end               
            end
            
            f1=1-sqrt(1-e1s);
            f2=1-sqrt(1-e2s);
            
            
            b1=a1.*(1-f1);
            b2=a2.*(1-f2);
            
            
            r1n=(a1.*a1.*b1).^(1/3);
            
            r2n=(a2.*a2.*b2).^(1/3);
            
            
            delta_r1=r1-r1n;
            delta_r2=r2-r2n;
            
            a1=a1+delta_r1;
            a2=a2+delta_r2;
            
            eps2=max(abs([delta_r1 delta_r2]));
            i2=i2+1;
            
            if (i2>10)
                break
            end
            
        end
        
        
        fh(2,:)=f2;
        fh(1,:)=f1;
        
    case 'RD'
        
        M2=4/3*pi*rho2*r2.^3;
        M1=4/3*pi*(rho1-rho2).*r1.^3;
        
        M=M1+M2;
        
        Ch1=0.4*(M1).*(r1.^2);
        Ch2=0.4*(M2).*(r2.^2);
        
        Ch=Ch1+Ch2;
        
        lambdah=Ch./(M.*r2.^2);
        
        Q=2.5*(1-1.5*lambdah);
        
        kf=(4-Q.^2)./(1+Q.^2);
        
        q=(r2.^3).*(Omega.^2)./(G.*M);
        
        w=q.*kf;
        
        fh(2,:)=3*w./(6+w);
        fh(1,:)=NaN;
        
        
    otherwise
        error('Unsupported order');
        
end
            
            