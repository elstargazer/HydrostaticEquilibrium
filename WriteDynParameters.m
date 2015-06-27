

filename = 'DynParameters.txt';

[a2,b2,c2] = fr2abc(rcorei,fp2_noneq,fq2_noneq);
[a2eq,b2eq,c2eq] = fr2abc(rcorei,fcorei{ModelIn},0);


in = fopen(filename,'w');

fprintf(in,'M = %23.16E [kg],\nV = %23.16E [m3],\nomega = %23.16E [rad/s],\na1 = %23.16E [m],\nb1 = %23.16E [m],\nc1 = %23.16E [m]\n',...
       M,V(ModelIn),W,a3*1000,b3*1000,c3*1000);
  
Ngrid=numel(rhoouteri{ModelIn});


Condition = ((rhoouteri{ModelIn} < 1500) & (rhocorei < 3500) & (rhoouteri{ModelIn} > 800));

fprintf(in,'a2noneq [m], b2noneq [m], c2noneq [m], a2eq [m], b2eq [m], c2eq [m], rho1 [kg/m3], rho2 [kg/m3]\n');

for i=1:Ngrid
    if (Condition(i))
        fprintf(in,'%23.16E, %23.16E, %23.16E, %23.16E, %23.16E, %23.16E, %23.16E, %23.16E\n',...
            a2(i),b2(i),c2(i),a2eq(i),b2eq(i),c2eq(i),rhoouteri{ModelIn}(i),rhocorei(i));
    end
end

fclose(in);