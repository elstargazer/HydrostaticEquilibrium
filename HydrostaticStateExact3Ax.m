function [fh,fval]=HydrostaticStateExact3Ax(r,T,rho,fp0,fq0)

A = [-1 1];
b = 0;

fp=0.001:0.02:0.9;
fq=0.0011:0.02:0.9;

[fp,fq]=meshgrid(fp,fq);
d2 = zeros(1,numel(fp));

progressbar(0);
parfor i=1:numel(fp)
    
    d2(i)=DeltaSquared3Ax([fp(i) fq(i)],r,T,rho);
%     progressbar(i/numel(fp));
end
progressbar(1);

d2=reshape(d2,size(fp));

fig1 = figure; hold on;
pcolor(fp,fq,log10(d2)); shading interp
caxis([4 10]);
contour(fp,fq,log10(d2),60,'Color','k');
shading interp
% Get coordinated of the minimum
mincoord=ginput(1);
close(fig1);

[fh,fval]=fminsearchcon(@(f) DeltaSquared3Ax(f,r,T,rho),...
    mincoord,[0 0],[1 1],A,b,[],...
    optimset('TolFun',1e-10,'MaxFunEvals',100000));



