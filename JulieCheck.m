fntsize=20;

[a,c]=rf2ac(router(1),fouteri{1});

R=MeanEllRadius(a,a,c);
Rv=(a.*a.*c).^(1/3);
Ra=(a+a+c)/3;

q=(a-R)./(c-R);
qv=(a-Rv)./(c-R);
qa=(a-Ra)./(c-Ra);

level_set2=0:250:6000;

figure; hold on;
set(gca,'FontSize',fntsize,'Position',[0.15 0.15 0.7 0.7]);
pcolor(rcorei/1000,rhocorei,real(q)); shading interp
cbar=colorbar('FontSize',fntsize);
ylabel(cbar,' (a-R_{mean})/(c-R_{mean}) ','FontSize',fntsize);
xlabel('Core size [km]','FontSize',fntsize);
ylabel('Core density [kg/m^3]','FontSize',fntsize);
xlim([0 470]);
ylim([2000 6000]);
caxis([-.54 -0.49])

figure; hold on;
set(gca,'FontSize',fntsize,'Position',[0.15 0.15 0.7 0.7]);
pcolor(rcorei/1000,rhocorei,real(qv)); shading interp
cbar=colorbar('FontSize',fntsize);
ylabel(cbar,' (a-R_{vol})/(c-R_{vol}) ','FontSize',fntsize);
xlabel('Core size [km]','FontSize',fntsize);
ylabel('Core density [kg/m^3]','FontSize',fntsize);
xlim([0 470]);
ylim([2000 6000]);
caxis([-.54 -0.49])

figure; hold on;
set(gca,'FontSize',fntsize);
pcolor(rcorei/1000,rhocorei,qa); shading interp
cbar=colorbar('FontSize',fntsize);
ylabel(cbar,' (a-R_{a})/(c-R_{a}) ','FontSize',fntsize);
xlabel('Core size [km]','FontSize',fntsize);
ylabel('Core density [kg/m^3]','FontSize',fntsize);

% contour(rcorei/1000,rhocorei,fouteri{1})

