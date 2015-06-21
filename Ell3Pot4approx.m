function U=Ell3Pot4approx(a,b,c,x,y,z,rho)

G=6.67e-11;

eq2=1-(b/a)^2;
ep2=1-(c/a)^2;


r=sqrt(x*x+y*y+z*z);

U=pi*G*rho*(2/3*(3*a*a-r*r)-2/15*(5*a*a-r*r+3*z*z)*ep2-2/15*(5*a*a-r*r+3*y*y)*eq2-...
    2/105*(14*a*a-4*r*r+12*z*z)*ep2*ep2+2/105*(14*a*a+r*r-3*x*x)*ep2*eq2-2/105*(14*a*a-4*r*r+12*y*y)*eq2*eq2);