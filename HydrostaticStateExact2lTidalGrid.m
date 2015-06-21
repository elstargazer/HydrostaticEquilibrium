function [fp1,fp2,fq1,fq2,fvali]=...
    HydrostaticStateExact2lTidalGrid(r1,r2,T,rho1,rho2)

% progressbar(0)

fhi=zeros(size(r2,1),4);
fvali=zeros(size(r2,1),1);

M=4;

progressbar(0)

parfor i=1:numel(rho2)
    if ((rho1(i)>0) && (rho2(i)>0) && (r2(i)>0))
        
        fh0=HydrostaticState2LayerTidalAn(r1,r2(i),...
            T,rho1(i),rho2(i),M);
        
        try     
            [fhi(i,:),fvali(i)]=...
                HydrostaticStateExact2lTidal(r1,r2(i),T,...
                rho1(i),rho2(i),fh0);
        catch
            fhi(i,:)=[NaN NaN NaN NaN];
            fvali(i)=NaN;
        end
    else
        fhi(i,:)=[NaN NaN NaN NaN];
        fvali(i)=NaN;
    end
%     progressbar(i/numel(rho2));
end

progressbar(1);

fp1=fhi(:,1);
fp2=fhi(:,2);
fq1=fhi(:,3);
fq2=fhi(:,4);

fp1=reshape(fp1,size(rho2));
fp2=reshape(fp2,size(rho2));
fq1=reshape(fq1,size(rho2));
fq2=reshape(fq2,size(rho2));

fvali=reshape(fvali,size(rho2));