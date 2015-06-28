function [fp, fq, fabp] = ell2f(ell)

fp = (ell(1)-ell(3))/ell(1);
fq = (ell(1)-ell(2))/ell(1);

ab = sqrt(ell(1)*ell(2));
fabp = (ab - ell(3))/ab;