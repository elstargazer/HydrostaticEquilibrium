ccc

%% Input parameters


r=6371000;
rho=5517;


T0=2;
T1=4;
Tstep=0.05;


% period
T=T1:-Tstep:T0;

% T=2.8;


%% Initial guess

fstep=0.005;

maxfq=0.99;
maxfp=0.99;


fqi=0.00001:fstep:maxfq;
fpi=0.00001:fstep:maxfp;

[fqi,fpi]=meshgrid(fpi,fqi);

%  progressbar(0);

vidObj = VideoWriter('DeltaSq_HomoEarth_HD.avi');
open(vidObj);

movie_fig=figure('Position',[1 1 1400 1400]); hold on;
set(gca,'FontSize',20);
xlim([0 1]);
ylim([0 1]);

xlabel('Polar flattening []','FontSize',20); 
ylabel('Equatorial flattening []','FontSize',20); 

box on;


for j=1:numel(T)

    clear d2;
    d2=zeros(1,numel(fpi));
    
    Tc=T(j);

    parfor i=1:numel(fpi)
         if (fpi(i)>fqi(i))
            d2(i)=DeltaSquared3Ax([fpi(i) fqi(i)],r,Tc,rho);       
         else
             d2(i)=NaN;
         end    
%          progressbar(i/numel(fpi));
    end

d2=reshape(d2,size(fpi));

pcolor(fpi,fqi,log10(d2));
contour(fpi,fqi,log10(d2),60,'Color','k'); shading interp;

% [in1,in2]=localMaximum(-log10(d2),[tstep tstep]);
% img=-(log10(d2)-max(log10(d2(:))));
% [cent]=FastPeakFind(img,0,[],0,2);
% for k=1:2:numel(cent)/2
%     plot(fqi(cent(k),cent(k+1)),fpi(cent(k),cent(k+1)),'ow','MarkerSize',10,'MarkerFaceColor','k');
% end


title(['T=' num2str(T(j)) ' [hr]']);
cbar=colorbar;
ylabel(cbar,'log_{10}(\Delta^{2}) [m^{4} sec ^{-4}]','FontSize',20)
caxis([6 16]);

drawnow;

currFrame = getframe(gcf);
writeVideo(vidObj,currFrame);

unplot(2);
% unplot(numel(cent)/2)

end

close(vidObj);

% progressbar(1);
% 
% d2=reshape(d2,size(fpi));
% 
% 
% figure; hold on;
% 
% pcolor(fpi,fqi,d2);
% contour(fpi,fqi,log10(d2),60,'Color','k');
% shading interp
% colorbar

% Get coordinated of the minimum

mincoord=ginput(1);

% Find initial equilibrium
[f0,fval]=HydrostaticStateExact3Ax(r,T1,rho,mincoord(1),mincoord(2))


plot(f0(1),f0(2),'*w');


%% Main loop

clear d2

progressbar(0);

tic
for i=1:numel(T)
       
    [fh(i,:),d2(i)]=HydrostaticStateExact3Ax(r,T(i),rho,f0(1),f0(2));
    f0=fh(i,:);
    
    progressbar(i/numel(T));
    
end
toc
% 
% progressbar(1);
% 
% 
% %% Plot
% % try
% % figure(figuref)
% % catch
% % end
% 
figure;
hold on;
%
plot(T,fh(:,1),'r-.','LineWidth',2);
plot(T,fh(:,2),'b-.','LineWidth',2);
% 
% % figure; hold on
% % data=load('SixOrderHydroList');
% % 
% % plot(data(3,:),data(1,:),'g');
% % plot(data(3,:),data(2,:),'r');
% 
% 
% figure;
% semilogy(T,d2,'b');
% grid on;
% 
