function [fp,fq,fabp,fp_err]=abc2f(abc,abc_err)

a=abc(1);
b=abc(2);
c=abc(3);

sa=abc_err(1);
sb=abc_err(2);
sc=abc_err(3);

fp = (a-c)./a;
fq = (a-b)./a;

fabp=(sqrt(a.*b)-c)./sqrt(a.*b);

fp_err = sqrt((c.*c.*sa.*sa + a.*a.*sc.*sc)./...
    (a.^4));