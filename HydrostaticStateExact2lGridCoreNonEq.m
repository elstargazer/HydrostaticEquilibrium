function [fp2,fq2,fval,exitflag]=...
    HydrostaticStateExact2lGridCoreNonEq(...
    router,rcorei,T,rhoouteri,rhocorei,fp1,fq1)

progressbar(0)

parfor i=1:numel(rhocorei)   

    if (rhoouteri(i)>0)

        
        [fhi(i,:),fval(i),exitflag(i)]=...
            HydrostaticStateExact2lCoreNonEq(...
            router,rcorei(i),T,rhoouteri(i),rhocorei(i),fp1,fq1);     
    else
        fhi(i,:) = [NaN NaN]; 
        fval(i) = NaN;
        exitflag(i) = NaN;
    end 
     % progressbar(i/numel(rhocorei));
end

progressbar(1);

fp2=fhi(:,1);
fq2=fhi(:,2);

fp2=reshape(fp2,size(rhocorei));
fq2=reshape(fq2,size(rhocorei));

fval=reshape(fval,size(rhocorei));
exitflag=reshape(exitflag,size(rhocorei));