function lmcosi=ReadSHNimmo(filename)

lmcosi_temp = load(filename);
CSness = lmcosi_temp(:,3);

C=lmcosi_temp(CSness == 1,4);
S=lmcosi_temp(CSness ==-1,4);

n=lmcosi_temp(CSness == 1,1);
m=lmcosi_temp(CSness == 1,2);


lmcosi = [0 0 1 0;
          1 0 0 0;
          1 1 0 0;
          n m C C*0];

k=1;
for n=2:lmcosi(end,1)
    for m=1:n
        
        ind = (n+1)*n/2+1+m;
        lmcosi(ind,4)=S(k);
        k=k+1;
    end 
end
