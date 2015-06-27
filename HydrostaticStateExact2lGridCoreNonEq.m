function [fp2,fq2,fval,exitflag]=...
    HydrostaticStateExact2lGridCoreNonEq(...
    router,rcorei,T,rhoouteri,rhocorei,fp1,fq1)

progressbar(0)

for i=1:numel(rhocorei)
    
    if (rhoouteri(i)>0)

        [fhi(i,:),fval(i),exitflag(i)]=...
            HydrostaticStateExact2lCoreNonEq(...
            router,rcorei(i),T,rhoouteri(i),rhocorei(i),fp1,fq1,0.1,0.05);
        
%         PlotPot2El(router,rcore,...
%             fp1,fq1,fhi(i,1),fhi(i,2),rhoouteri(i),rhocorei(i),W);
        
        
%         fouter0=fhi(i,1);
%         fcore0=fhi(i,2);      
    else
        fhi(i,:)=[NaN NaN]; 
        fval(i)=NaN;
        exitflag(i)=NaN;
    end 
%      progressbar(i/numel(rhocorei));
%      i/numel(rhocorei)
end

progressbar(1);

fp2=fhi(:,1);
fq2=fhi(:,2);

fp2=reshape(fp2,size(rhocorei));
fq2=reshape(fq2,size(rhocorei));

fval=reshape(fval,size(rhocorei));
exitflag=reshape(exitflag,size(rhocorei));