ccc;

M=4;

r1=252000;
r2=100;
rho1=1601;
rho2=2000;

Rref=252000;

Tn=24:1:130;

fntsize=20;
fntsize_sm=18;
imsize=[54 38];

%% Compute flattenings

progressbar(0);

parfor i=1:numel(Tn)
    
    [fh(i,:)]=HydrostaticState2LayerTidalAn...
        (r1,r2,Tn(i),rho1,rho2,M);
    
    [fhn(i,:),~]=HydrostaticStateExact2lTidal...
        (r1,r2,Tn(i),rho1,rho2,fh(i,:));
    
%     progressbar(i/numel(Tn));
end

progressbar(1);


fp1t=fh(:,1);
fp2t=fh(:,2);
fq1t=fh(:,3);
fq2t=fh(:,4);

fp1tn=fhn(:,1);
fp2tn=fhn(:,2);
fq1tn=fhn(:,3);
fq2tn=fhn(:,4);

%% Gravity harmonics J2 and C22

[J2,C22]=RadFlat2J2Tri(r1,r2,...
    fp1t,fp2t,fq1t,fq2t,rho1,rho2,Rref);

[J2n,C22n]=RadFlat2J2Tri(r1,r2,...
    fp1tn,fp2tn,fq1tn,fq2tn,rho1,rho2,Rref);

fg=J2./C22;
fgn=J2n./C22n;

fst=fp1t./(fp1t-fq1t);
fstn=fp1tn./(fp1tn-fq1tn);


%% Flattening plot
figure; hold on;
set(gca,'FontSize',fntsize);

plot(Tn,fp1t,'r-');
plot(Tn,fp2t,'g-');
plot(Tn,fq1t,'b-');
plot(Tn,fq2t,'m-');

plot(Tn,fp1tn,'r--');
plot(Tn,fp2tn,'g--');
plot(Tn,fq1tn,'b--');
plot(Tn,fq2tn,'m--');

xlabel('Period [h]','FontSize',fntsize);
ylabel('Flattening []','FontSize',fntsize);

%% 10/3 plot
figure; hold on;
set(gca,'FontSize',fntsize);

plot(Tn,J2./C22,'r-');
plot(Tn,J2n./C22n,'g-');

xlabel('Period [h]','FontSize',fntsize);
ylabel('J2/C22 []','FontSize',fntsize);

legend({'An','Num'},'FontSize',fntsize);


%% 1/4 plot
figure; hold on;
set(gca,'FontSize',fntsize);

plot(Tn,fst,'r-');
plot(Tn,fstn,'g-');

xlabel('Period [h]','FontSize',fntsize);
ylabel('J2/C22 []','FontSize',fntsize);

legend({'An','Num'},'FontSize',fntsize);
