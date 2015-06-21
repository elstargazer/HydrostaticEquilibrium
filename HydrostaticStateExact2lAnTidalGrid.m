function [fp1,fp2,fq1,fq2]=...
    HydrostaticStateExact2lAnTidalGrid(r1,r2,T,rho1,rho2,M)

Npts=numel(r2);
progressbar(0)
fh=zeros(4,Npts);
s=size(r2);

progressbar(0);
for i=1:Npts
    fh(:,i)=HydrostaticState2LayerTidalAn(r1,r2(i),...
        T,rho1(i),rho2(i),M);     
    progressbar(i/Npts); 
end
progressbar(1);

fp1=fh(1,:);
fp2=fh(2,:);
fq1=fh(3,:);
fq2=fh(4,:);

fp1=reshape(fp1,s);
fp2=reshape(fp2,s);
fq1=reshape(fq1,s);
fq2=reshape(fq2,s);
