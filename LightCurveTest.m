ccc

a=960.000;
b=770.000;
c=495.000;

r=715;

fq=(a-b)/a;
fp=(a-c)/a;

a=r/(((fp-1)*(fq-1))^(1/3));
b=(r - fq*r)/(((fp-1)*(fq-1))^(1/3));
c=(r - fp*r)/(((fp-1)*(fq-1))^(1/3));

area_ab_mag=-2.5*log10(pi*a*b);

area_cb_mag=-2.5*log10(pi*c*b)-area_ab_mag;

area_ca_mag=-2.5*log10(pi*c*a)-area_ab_mag;


alpha=4;


% [area]=EllipsoidCrossArea(r,fq,fp,fi,lambda)

[area_obs]=EllipsoidLightCurve(r,fq,fp,alpha);
area_obs_mag=-2.5*log10(area_obs)-area_ab_mag;

figure; hold on;

phase=linspace(0,1,numel(area_obs));

xlabel('Phase []','FontSize',12);
ylabel('Projected area -2.5 log_{10}(km^{2})','FontSize',12);

set(gca,'YDir','reverse');
box on;

plot(phase,0.*ones(size(area_obs)),'r-');
plot(phase,area_cb_mag.*ones(size(area_obs)),'g-');
plot(phase,area_ca_mag.*ones(size(area_obs)),'b-');
plot(phase,area_obs_mag,'-k','LineWidth',3);


legend({'a-b area','c-b area','c-a area','Light curve'},'FontSize',12);

fpi=fp/1.25:0.01:fp*1.25;
fqi=fq/1.25:0.01:fq*1.25;

for i=1:numel(fpi)
    [area]=EllipsoidLightCurve(r,fq,fpi(i),alpha);
    area_mag=-2.5*log10(area)-area_ab_mag;
    plot(phase,area_mag,'--m','LineWidth',1);
end


for i=1:numel(fqi)
    [area]=EllipsoidLightCurve(r,fqi(i),fp,alpha);
    area_mag=-2.5*log10(area)-area_ab_mag;
    plot(phase,area_mag,'Color',[0.8 0.7 0],'LineWidth',1);
end



fpi=0:0.008:0.999;
fqi=0:0.008:0.999;

[fpi,fqi]=meshgrid(fpi,fqi);

s=size(fqi);

fpi=fpi(:);
fqi=fqi(:);

parfor i=1:numel(fpi)
    
    [areai]=EllipsoidLightCurve(r,fqi(i),fpi(i),alpha);
    areai_mag=-2.5*log10(abs(areai))-area_ab_mag;
    
    v2(i)=sum(((areai_mag-max(areai_mag))-(area_obs_mag-max(area_obs_mag))).^2);  
end


fpi=reshape(fpi,s);
fqi=reshape(fqi,s);
v2=reshape(v2,s);

figure; hold on;

pcolor(fpi,fqi,log10(v2)); shading interp;

xlabel('Polar flattening []','FontSize',12);
ylabel('Equatorial flattening []','FontSize',12);

plot(fp/1.02,fq,'k*','MarkerSize',10);

cbar=colorbar;

ylabel(cbar,'\chi^{2}','FontSize',12);
box on;


    






























