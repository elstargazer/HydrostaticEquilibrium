ccc

G=6.67e-11;

GM=2.54e9;
T=23:-0.5:9;
M=GM/G;

GMp=(37.931e6)*1e9;

% T=0.942*86400
aorbit=(((T*3600).^2)*GMp/4/pi/pi).^(1/3);

a=207.98;
b=196.22;
c=191.40;

fp_obs=(a-c)/a;
fq_obs=(a-b)/a;

fs_obs=(a-c)/(b-c);

r1=198000;
V=4/3*pi*r1.^3;
rhomean=M/V;

fstep=0.0025;
maxfq=0.6;
maxfp=0.6;

fqi=0.00001:fstep:maxfq;
fpi=0.00001:fstep:maxfp;

[fqi,fpi]=meshgrid(fpi,fqi);

figure('Position',[1 1 1400 1400]); hold on;
set(gca,'FontSize',20);

writerObj = VideoWriter('HomoTidalD2.avi');
open(writerObj);

for j=1:numel(T)
    
    Tcurrent=T(j);
    parfor i=1:numel(fpi)
        if (fpi(i)>fqi(i))
%             d2(i)=DeltaSquared3Ax([fpi(i) fqi(i)],r1,Tcurrent,rhomean);
            d2(i)=DeltaSquaredTidal([fpi(i) fqi(i)],r1,Tcurrent,rhomean);     
        else
            d2(i)=NaN;
        end
        % progressbar(i/numel(fpi));
    end
    
    Omega=2*pi./(Tcurrent*3600);
    Lambda2=(Omega.^2)./(pi*G*rhomean);
    
    fsinv = 0.25 - 1485/896 * Lambda2;
    fs=1./fsinv;
    
    d2=reshape(d2,size(fpi));
    pcolor(fpi,fqi,log10(d2)); shading interp;
    levels=-40:.25:40;
    contour(fpi,fqi,log10(d2),levels,'Color','k'); 
    contour(fpi, fqi, fpi./(fpi-fqi),[4 4],'Color','w','LineWidth',2);
    contour(fpi, fqi, fpi./(fpi-fqi),[fs fs],'Color','g','LineWidth',2);
    plot(fp_obs,fq_obs,'om','MarkerSize',9);
    
    cbar=colorbar;
    ylabel(cbar,'log_{10}(\Delta^{2}) [m^{4} sec ^{-4}]','FontSize',20)
    
    xlim([0 maxfq]);
    ylim([0 maxfp]);
    
    xlabel('Polar flattening []','FontSize',20);
    ylabel('Equatorial flattening []','FontSize',20);
    caxis([min(log10(d2(:))) max(log10(d2(:)))]);
    title(['T = ' num2str(Tcurrent) ' [h]; a = ' num2str(aorbit(j)/1000) ...
    ' [km]'],'FontSize',20);
    box on
    drawnow
    frame = getframe;
    writeVideo(writerObj,frame);
    
end

close(writerObj);