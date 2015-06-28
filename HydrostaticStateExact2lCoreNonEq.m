function [fh,fval,exitflag]=HydrostaticStateExact2lCoreNonEq(r1,r2,T,rho1,rho2,fp1,fq1,fp2_0,fq2_0)

A=[-1 1];
b=0;

% [fh,fval]=fminsearch(@(f) DeltaSquared2lCoreNonEq(f,r1,r2,T,rho1,rho2,fp1,fq1),[f10 f20],...
%     optimset('TolFun',1e-10));

% 
<<<<<<< HEAD
fp2i=linspace(0.001,0.5,20);
fq2i=linspace(0.001,0.5,20);
=======
fp2i=linspace(0.0001,0.9,30);
fq2i=linspace(0.00011,0.9,30);

>>>>>>> 3a0ebf8cff0cd8c2f9067f3e09b40bc6954ce854

[fq2i,fp2i]=meshgrid(fq2i,fp2i);
d2=nan(size(fq2i));

for i=1:numel(fq2i)
<<<<<<< HEAD
        d2(i)=DeltaSquared2lCoreNonEq(...
            [fp2i(i) fq2i(i)],...
            r1,r2,T,rho1,rho2,...
            fp1,fq1); 
=======
    
    d2(i)=DeltaSquared2lCoreNonEq(...
        [fp2i(i) fq2i(i)],...
        r1,r2,T,rho1,rho2,...
        fp1,fq1);
    
>>>>>>> 3a0ebf8cff0cd8c2f9067f3e09b40bc6954ce854
end

index=find(d2==min(d2(:)));

fq2init=fq2i(index);
fp2init=fp2i(index);

% d2 = reshape(d2,size(fp2i));

% figure; hold on;
<<<<<<< HEAD
% surf(fp2i,fq2i,d2-min(d2(:))); shading interp
% contour(fp2i,fq2i,log10(d2),100,'Color','k'); 

%  [fh,fval,exitflag,~]=fminsearchcon(@(f) DeltaSquared2lCoreNonEq(f,r1,r2,T,rho1,rho2,fp1,fq1),...
%      [fp2init fq2init],[0 0],[1 1],A,b,[],...
%      optimset('MaxFunEvals',1e5,'TolX',1e-4));
 
eps = 1e-6;

  [fh,fval,exitflag,~]=fminsearchcon(@(f) DeltaSquared2lCoreNonEq(f,r1,r2,T,rho1,rho2,fp1,fq1),...
     [fp2init fq2init],[0 0],[1 1]-eps,[],[],[],...
     optimset('MaxFunEvals',1e5,'TolX',1e-9,'TolFun',1e-7));
 
 
=======
% pcolor(fp2i,fq2i,log10(d2)); shading interp
% contour(fp2i,fq2i,log10(d2),50,'Color','k'); 
% colorbar;

 [fh,fval,exitflag,~]=fminsearchcon(@(f) ...
     DeltaSquared2lCoreNonEq(f,r1,r2,T,rho1,rho2,fp1,fq1),...
     [fp2init fq2init],[0 0],[1 1],[],[],[],...
     optimset('MaxFunEvals',1e5,'TolX',1e-4));
 
 fval
>>>>>>> 3a0ebf8cff0cd8c2f9067f3e09b40bc6954ce854
