function [sr, sr_udg] = sourcema2(p,udg,param,time)

R = param(1);
x = p(:,1);
y = p(:,2);
sr =  -(R^2)./(R^2 - x.^2 - y.^2).^2;

[ng,nc] = size(udg);
nch = 1;
sr_udg = zeros(ng,nch,nc); 
