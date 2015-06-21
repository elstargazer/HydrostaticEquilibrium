function [area]=EllipsoidLightCurve(r,fq,fp,alpha)

[latc,lonc] = scircle1(89.999,0,90-alpha);


% figure;
% 
% plot(lonc,latc);


[area]=EllipsoidCrossArea(r,fq,fp,latc/180*pi,lonc/180*pi);