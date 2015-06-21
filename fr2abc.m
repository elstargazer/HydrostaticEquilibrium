function [a,b,c]=fr2abc(r,fp,fq)

a=r./(((fp-1).*(fq-1)).^(1/3));
b=(r - fq.*r)./(((fp-1).*(fq-1)).^(1/3));
c=(r - fp.*r)./(((fp-1).*(fq-1)).^(1/3));