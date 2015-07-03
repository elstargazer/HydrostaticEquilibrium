ccc
rho=1000;
a=500000;
b=450000;
c=400000;

mu=muEllipsoid(a,b,c,rho);

Nmax = 2:2:60;

tic
parfor j=1:numel(Nmax)
    
    N=Nmax(j);
    N_trunc=Nmax(j);
    Rref=500000;
    
    [lmcosi]=SHTriaxialEllipsoid(a,b,c,N,N_trunc,Rref);
    
    TRr=IcosahedronMesh;
    TR_2r=SubdivideSphericalMesh(TRr,2);
    [lambda,fi,~]=cart2sph(TR_2r.X(:,1),TR_2r.X(:,2),TR_2r.X(:,3));
    [x,y,z]=TriEllRadVec(fi,lambda,a,b,c,'xyz');
    
    Uexact=zeros(numel(x),1);
    for i=1:numel(x)
        Uexact(i)=Ell3Pot(a,b,c,x(i),y(i),z(i),rho);
    end
    
    Ush = GravityPotential(mu,Rref,lmcosi,x,y,z);    
    n2(j) = sqrt(sum((Uexact-Ush).^2)/numel(Ush));
    
end
toc

figure; hold on;
set(gca,'YScale','log');
plot(Nmax,n2,'ok','MarkerSize',6);
xlabel('Max SH degree','FontSize',12);
ylabel('Norm of diff','FontSize',12);








