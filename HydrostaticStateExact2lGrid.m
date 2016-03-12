function [fcorei,fouteri]=HydrostaticStateExact2lGrid...
    (router,rcorei,T,rhoouteri,rhocorei)

progressbar(0)

parfor i=1:numel(rhocorei)   
    if (rhoouteri(i)>0 && (rhoouteri(i)~=Inf))
        
        [rcorei(i) rhoouteri(i),rhocorei(i)]
        [fhi(i,:),~]=HydrostaticStateExact2l(router,rcorei(i),T,...
            rhoouteri(i),rhocorei(i),0.1,0.1);    
    else
        fhi(i,:)=[NaN NaN];       
    end  
%     progressbar(i/numel(rhocorei));
end

progressbar(1);

fcorei=fhi(:,2);
fouteri=fhi(:,1);

fcorei=reshape(fcorei,size(rhocorei));
fouteri=reshape(fouteri,size(rhocorei));