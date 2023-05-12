function [sr, sr_udg] = sourcema4(p,udg,param,time)

x = p(:,1);
y = p(:,2);
sr = -exp(x.^2 + y.^2).*(x.^2 + y.^2 + 1);

v11 = p(:,3);
v12 = p(:,4);
v21 = p(:,5);
v22 = p(:,6);

sr = sr - v11.*v22 + v12.*v21;

[ng,nc] = size(udg);
nch = 1;
sr_udg = zeros(ng,nch,nc); 
