% ccc

M=4.006e21;

T=3.915341;

a_guess=960000;
b_guess=770000;
c_guess=495000;

r1=(a_guess*b_guess*c_guess)^(1/3);
V=4/3*pi*r1^3;

rhomean=M/V;
 
 %% Two layer model
 
rhocore=rhomean:50:5000;
rcore=100:5000:r1;
router=r1;

[rhocorei,rcorei]=meshgrid(rhocore,rcore);

Vcore=4/3*pi*rcorei.^3;
Mcore=rhocorei.*Vcore;

Vmantle=V-Vcore;
Mmantle=M-Mcore;

rhomantle=Mmantle./Vmantle;

rhoouter=MantleDensity(M,rcorei,rhocorei,router);

l2m_fig=figure; hold on;

pcolor(rcorei/1000,rhocorei,rhomantle);

levels=[0:500:3000 3100:100:rhomean];
[C,h]=contour(rcorei/1000,rhocorei,rhomantle,levels,'Color','k');
text_handle = clabel(C,h);
set(text_handle,'BackgroundColor',[1 1 .6],...
    'Edgecolor',[.7 .7 .7])

cbar=colorbar;
ylabel(cbar,'Outer density [kg/m^{3}]','FontSize',12);
box on;

xlabel('Core radius [km]','FontSize',12);
ylabel('Core density [kg/m^{3}]','FontSize',12);

shading interp

xlim([0 r1/1000]);
ylim([rhomean 5000]);
caxis([500 rhomean]);


%% Hydrostatic equilibrium for homogeneous body

f_p_guess=(a_guess-c_guess)/a_guess;
f_q_guess=(a_guess-b_guess)/a_guess;

[f_p_guess f_q_guess]

fstep=0.005;

maxfq=0.99;
maxfp=0.99;

fqi=0.00001:fstep:maxfq;
fpi=0.00001:fstep:maxfp;

[fqi,fpi]=meshgrid(fpi,fqi);


parfor i=1:numel(fpi)
      if (fpi(i)>fqi(i))
         d2(i)=DeltaSquared3Ax([fpi(i) fqi(i)],r1,T,rhomean);       
      else
        d2(i)=NaN;
      end    
%          progressbar(i/numel(fpi));
end

d2=reshape(d2,size(fpi));

figure; hold on;
set(gca,'FontSize',20);


pcolor(fpi,fqi,log10(d2));
contour(fpi,fqi,log10(d2),60,'Color','k'); shading interp;

cbar=colorbar;
ylabel(cbar,'log_{10}(\Delta^{2}) [m^{4} sec ^{-4}]','FontSize',20)

xlim([0 1]);
ylim([0 1]);

xlabel('Polar flattening []','FontSize',20); 
ylabel('Equatorial flattening []','FontSize',20); 
caxis([6 16]);


f_homo_guess=ginput(1);

tic
[f0,fval]=HydrostaticStateExact3Ax(r1,T,rhomean,f_homo_guess(1),f_homo_guess(2))
toc

%% Hydrostatic equilibrium 2 layer

delta_rad=-10000;
delta_rho=10;

rcore_test=r1+delta_rad;
rhocore_test=rhomean+delta_rho;
rhoouter_test=MantleDensity(M,rcore_test,rhocore_test,r1);

figure(l2m_fig);

plot(rcore_test/1000,rhocore_test,'w*','MarkerSize',10);

f_test=[f0(2) f0(1) f0(2)-0.001 f0(1)-0.001];

tic
[f_initial,fval]=HydrostaticStateExact3Ax2l(r1,rcore_test,T,rhoouter_test,rhocore_test,f_test)
toc


% making grid
rad_step=7500;
rho_step=50;


rcore2=r1+delta_rad:-rad_step:100;
rhocore2=rhomean+delta_rho:rho_step:5000;

[rcore2i,rhocore2i]=meshgrid(rcore2,rhocore2);

figure(l2m_fig);

plot(rcore2i/1000,rhocore2i,'w*','MarkerSize',10);

fp1=zeros(size(rcore2i));
fp2=zeros(size(rcore2i));
fq1=zeros(size(rcore2i));
fq2=zeros(size(rcore2i));
fvali=zeros(size(rcore2i));

f_guess=f_initial;

tic

for i=1:size(rcore2i,1)
    
%     if (i>1)
%         f_guess=f_guess_fr;
%     end
    
    parfor j=1:size(rcore2i,2)
        
        f_guess=f_test;
        
            
        
        rhoouter2=MantleDensity(M,rcore2i(i,j),rhocore2i(i,j),r1);
        
        if (rhoouter2>0)
            tic
            [fh,fval]=HydrostaticStateExact3Ax2l(r1,rcore2i(i,j),T,...
                rhoouter2,rhocore2i(i,j),f_guess);
            toc
        else
            fh=nan(4,1)
            fval=NaN;
        end
        
        fq1(i,j)=fh(1);
        fp1(i,j)=fh(2);
        fq2(i,j)=fh(3);
        fp2(i,j)=fh(4);
        
        fvali(i,j)=fval;
        
        plot(rcore2i(i,j)/1000,rhocore2i(i,j),'m*','MarkerSize',10);
        drawnow
        
        f_guess=f_initial;
        
%         if (j==1)
%             f_guess_fr=fh;
%         end
     
    end
    
     i/size(rcore2i,1)*100
    
end

toc


% Plot residual
figure; hold on;

pcolor(rcore2i/1000,rhocore2i,log10(fvali)); shading flat;

xlabel('Core radius [km]','FontSize',12);
ylabel('Core density [kg/m^{3}]','FontSize',12);


fq1(log10(fvali)>-2)=NaN;
fp1(log10(fvali)>-2)=NaN;

rhoouter2i=MantleDensity(M,rcore2i,rhocore2i,r1);




% plotting fq1

levels=0.0:0.05:1;

figure; hold on;

subplot(2,2,1); hold on;
set(gca,'FontSize',15);


pcolor(rcore2i/1000,rhocore2i,fq1); shading flat;
[C,h]=contour(rcore2i/1000,rhocore2i,fq1,levels,'Color','k'); shading flat;
text_handle = clabel(C,h);
set(text_handle,'BackgroundColor',[1 1 .6],...
    'Edgecolor',[.7 .7 .7],'FontSize',12)
xlabel('Core radius [km]','FontSize',15);
ylabel('Core density [kg/m^{3}]','FontSize',15);
shading interp
colorbar
title('Outer equatorial flattening []','FontSize',15);
box on;

[C,h]=contour(rcore2i/1000,rhocore2i,rhoouter2i./rhocore2i,-10:0.2:10,'Color','r');
text_handle = clabel(C,h);
set(text_handle,'BackgroundColor',[1 1 .6],...
    'Edgecolor',[.7 0.0 .7],'FontSize',12)


caxis([0 1])



% plotting fp1

levels=0:0.025:1;

% figure; hold on;

subplot(2,2,2);hold on;
set(gca,'FontSize',15);

rhoouter2i=MantleDensity(M,rcore2i,rhocore2i,r1);



pcolor(rcore2i/1000,rhocore2i,fp1); shading flat;
[C,h]=contour(rcore2i/1000,rhocore2i,fp1,levels,'Color','k'); shading flat;
text_handle = clabel(C,h);
set(text_handle,'BackgroundColor',[1 1 .6],...
    'Edgecolor',[.7 .7 .7],'FontSize',12)
xlabel('Core radius [km]','FontSize',15);
ylabel('Core density [kg/m^{3}]','FontSize',15);
shading interp
colorbar

title('Outer polar flattening []','FontSize',15);
box on;

% plotting fq1

levels=0.0:0.05:1;

% figure; hold on;

subplot(2,2,3);hold on;
set(gca,'FontSize',15);


pcolor(rcore2i/1000,rhocore2i,fq2); shading flat;
[C,h]=contour(rcore2i/1000,rhocore2i,fq2,levels,'Color','k'); shading flat;
text_handle = clabel(C,h);
set(text_handle,'BackgroundColor',[1 1 .6],...
    'Edgecolor',[.7 .7 .7],'FontSize',12)
xlabel('Core radius [km]','FontSize',15);
ylabel('Core density [kg/m^{3}]','FontSize',15);
shading interp
colorbar
title('Inner equatorial flattening []','FontSize',15);
box on;


% plotting fp1

levels=0:0.025:1;

% figure; hold on;

subplot(2,2,4);hold on;
set(gca,'FontSize',15);



pcolor(rcore2i/1000,rhocore2i,fp2); shading flat;
[C,h]=contour(rcore2i/1000,rhocore2i,fp2,levels,'Color','k'); shading flat;
text_handle = clabel(C,h);
set(text_handle,'BackgroundColor',[1 1 .6],...
    'Edgecolor',[.7 .7 .7],'FontSize',12)
xlabel('Core radius [km]','FontSize',15);
ylabel('Core density [kg/m^{3}]','FontSize',15);
shading interp
colorbar

title('Inner polar flattening []','FontSize',15);
box on;

% figure; hold on;
% pcolor(rcore2i/1000,rhocore2i,fp1); shading flat;
% [C,h]=contour(rcore2i/1000,rhocore2i,fq1./fp1,levels,'Color','k'); shading flat;
% text_handle = clabel(C,h);
% set(text_handle,'BackgroundColor',[1 1 .6],...
%     'Edgecolor',[.7 .7 .7])
% xlabel('Core radius [km]','FontSize',12);
% ylabel('Core density [kg/m^{3}]','FontSize',12);
% shading interp


%% Differention

% fp_obs=f_p_guess;
% fq_obs=f_q_guess;

coords=ginput(1)

r2_guess=coords(1)*1000;
rho2_guess=coords(2);
r1_guess=r1;
fp1_guess=griddata(rcore2i,rhocore2i,fp1,r2_guess,rho2_guess,'cubic');
fq1_guess=griddata(rcore2i,rhocore2i,fq1,r2_guess,rho2_guess,'cubic');
fp2_guess=griddata(rcore2i,rhocore2i,fp2,r2_guess,rho2_guess,'cubic');
fq2_guess=griddata(rcore2i,rhocore2i,fq2,r2_guess,rho2_guess,'cubic');
rho1_guess=MantleDensity(M,r2_guess,rho2_guess,r1);




x_guess=[rho1_guess rho2_guess fq2_guess fp2_guess];

x_guess(1)
x_guess(2)
x_guess(3)
x_guess(4)


[xh,fval,exitflag,output]=HydrostaticStateExact3Ax2lmod(fq1_guess,fp1_guess,r1,r2_guess,T,x_guess);

xh'
fval'


[fh,fval]=HydrostaticStateExact3Ax2l(r1,r2_guess,T,rho1_guess, rho2_guess,...
    [fq1_guess fp1_guess fq2_guess fp2_guess]);

fh'
fval'

xh(1)
xh(2)
xh(3)
xh(4)

%% Grid differention

fp_obs=f0(1)/1.02;
fq_obs=f0(2);


r1i=r1-300000:40000:r1+300000;
r2i=20000:40000:max(r1i);


[r1i,r2i]=meshgrid(r1i,r2i);

s=size(r1i);

r1i=r1i(:);
r2i=r2i(:);

rho1i=zeros(size(r1i));
rho2i=zeros(size(r2i));
fp2i=zeros(size(r2i));
fq2i=zeros(size(r2i));
fvali=zeros(size(r2i));

numel(fvali)




parfor i=1:numel(r1i)
    if (r1i(i)>r2i(i))
        [xh,fval,exitflag,output]=HydrostaticStateExact3Ax2lmod(fq_obs,fp_obs,r1i(i),r2i(i),T,x_guess);
        
        rho1i(i)=xh(1);
        rho2i(i)=xh(2);
        fq2i(i)=xh(3);
        fp2i(i)=xh(4);
        fvali(i)=fval;
        fval
        
    else      
        rho1i(i)=NaN;
        rho2i(i)=NaN;
        fq2i(i)=NaN;
        fp2i(i)=NaN;
    end
%     i
end


rho1i=reshape(rho1i,s);
rho2i=reshape(rho2i,s);
fq2i=reshape(fq2i,s);
fp2i=reshape(fp2i,s);
fvali=reshape(fvali,s);
r1i=reshape(r1i,s);
r2i=reshape(r2i,s);


Mi=4/3*pi*(r1i.^3).*(rho1i)+4/3*pi*(r2i.^3).*(rho2i-rho1i);
Vi=4/3*pi*(r1i.^3);

rhomeani=Mi./Vi;


% Mass constraint plot

figure; hold on;

levels=0:0.1:10;
fvali_plot=pcolor(r1i/1000,rhomeani,log10(fvali)); shading flat;
colorbar;


figure; hold on;
% pcolor(r1i/1000,rhomeani,Mi);

pcolor(r1i/1000,rhomeani,log10(fvali)); shading flat;
colorbar


levels=[1 1];
[C,h]=contour(r1i/1000,rhomeani,Mi/M,levels,'Color','k','LineWidth',4);
text_handle = clabel(C,h);
set(text_handle,'BackgroundColor',[1 1 .6],...
    'Edgecolor',[.7 .7 .7])

xlabel('Radius [km]','FontSize',12);
ylabel('Mean density [kg/m^{3}]','FontSize',12);

xlim([min(r1i(:)) max(r1i(:))]/1000);
ylim([0 rhomean]);

%

figure; hold on;

levels=0:0.1:10;
fvali_plot=pcolor(r1i/1000,r2i/1000,log10(fvali)); shading flat;
cbar=colorbar;
ylabel(cbar,'Residual potential squared [m^2 sec^{-4}]','FontSize',12)
% 
% [C,h]=contour(r1i/1000,r2i/1000,Mi/M,'Color','k');
% text_handle = clabel(C,h);
% set(text_handle,'BackgroundColor',[1 1 .6],...
%     'Edgecolor',[.7 .7 .7])

xlabel('Radius outer [km]','FontSize',12);
ylabel('Radius inner [km]','FontSize',12);
box on;
        
%

figure; hold on;

% levels=0.0:0.1:1;
fvali_plot=pcolor(r1i/1000,r2i/1000,log10(fvali)); shading flat;
cbar=colorbar;
ylabel(cbar,'Residual potential squared [m^2 sec^{-4}]','FontSize',12)

[C,h]=contour(r1i/1000,r2i/1000,rhomeani,'Color','k');
text_handle = clabel(C,h);
set(text_handle,'BackgroundColor',[1 1 .6],...
    'Edgecolor',[.7 .7 .7])


[C,h]=contour(r1i/1000,r2i/1000,Mi/M,[1 1],'Color','k','LineWidth',4);

caxis([-6 4]);

xlabel('Radius outer [km]','FontSize',12);
ylabel('Radius inner [km]','FontSize',12);
        






        
        





















