ccc

f10=0.1;
f20=0.1;

rho1=4000;
rho2=2200;


r1 = 110000;
r2 = 476200;

T=9;

GM=62.8981e9;
G=6.67e-11;
M=GM/G;

% Order=2

HydrostaticState2LayerAn(r1,r2,T,rho1,rho2,f10,f20,1)

HydrostaticState2LayerAn(r1,r2,T,rho1,rho2,f10,f20,2)

HydrostaticState2LayerAn(r1,r2,T,rho1,rho2,f10,f20,6)

HydrostaticState2LayerAn(r1,r2,T,rho1,rho2,f10,f20,'RD')



T=8:0.1:10;

for i=1:numel(T)

    fh1(:,i)=HydrostaticState2LayerAn(r1,r2,T(i),rho1,rho2,f10,f20,1);
    fh2(:,i)=HydrostaticState2LayerAn(r1,r2,T(i),rho1,rho2,f10,f20,2);
    fh6(:,i)=HydrostaticState2LayerAn(r1,r2,T(i),rho1,rho2,f10,f20,6);

end


figure; hold on;

plot(T,fh1,'-r');
plot(T,fh2,'-g');
plot(T,fh6,'-b');


xlabel('Rotation period [h]','FontSize',12);
ylabel('Polar flattening []','FontSize',12);


min_r1=0;
max_r1=476000;

min_rho1=2200;
max_rho1=6000;

fobs=0.0669;
fobs_sigma=0.0055;
T=9.07


i=1;

Np=10000;

progressbar(0);

while (i<Np)

    r1=rand(1)*(max_r1-min_r1)+min_r1;
    rho1=rand(1)*(max_rho1-min_rho1)+min_rho1;

    rho2=-(3*M-4*pi*(r1.^3).*rho1)./(4*pi*(r1.^3)-4*pi*(r2^3));

    if (rho2<rho1)
        fh=HydrostaticState2LayerAn(r1,r2,T,rho1,rho2,f10,f20,1);

        if ((fh(1)>0) && (fh(2)>0))
            
            fhi(:,i)=fh;
        
            r1i(i)=r1;
            rho1i(i)=rho1;

            i=i+1;
        end

      
    end

progressbar(i/Np);


end

progressbar(1);


chi2=((fhi(2,:)-fobs).^2)./(fobs_sigma.^2);

figure; hold on;

xlim([min_r1 max_r1]/1000);
ylim([min_rho1 max_rho1]);

xlabel('Core radius [km]','FontSize',12);
ylabel('Core density [kg/m^{3}]','FontSize',12);

scatter(r1i/1000,rho1i,15,log10(chi2),'filled');

colorbar

chi2















