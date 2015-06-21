function k=GetK(a,b,c,x,y,z)

r2=x.*x+y.*y+z.*z;

k=r2-a*a;

Tol=1e-14;

d=Tol+1;

while (d>Tol)
    
    
    k0=k;
    k=x.*x./(1+a.*a./k)+y.*y./(1+b.*b./k)+z.*z./(1+c.*c./k);
    
    d=abs((k-k0)/k);   
    
    
end

% k