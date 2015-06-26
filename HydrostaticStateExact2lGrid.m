function [fcorei,fouteri]=HydrostaticStateExact2lGrid(router,rcorei,T,rhoouteri,rhocorei)

progressbar(0)

for i=1:numel(rhocorei)
    
    if (rhoouteri(i)>0)

        [fhi(i,:),fval]=HydrostaticStateExact2l(router,rcorei(i),T,...
            rhoouteri(i),rhocorei(i),0.1,0.1);
        %[fhi(i,:)]=HydrostaticState2LayerAn(rcorei(i),router,T,rhocorei(i),rhoouteri(i),0.1,0.1,1);
        
        fouter0=fhi(i,1);
        fcore0=fhi(i,2);      
    else
        fhi(i,:)=[NaN NaN];       
    end
    
    progressbar(i/numel(rhocorei));
end

progressbar(1);

fcorei=fhi(:,2);
fouteri=fhi(:,1);

fcorei=reshape(fcorei,size(rhocorei));
fouteri=reshape(fouteri,size(rhocorei));