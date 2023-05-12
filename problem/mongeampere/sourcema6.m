function [sr, sr_udg] = sourcema6(p,udg,param,time)

theta = param(1);
a1 = param(2);
a2 = param(3);
a = param(4);

qx = udg(:,2);
qy = udg(:,3);
r = sqrt((qx-1.5).^2 + qy.^2);
rho = 1 + a1*sech(a2*(r.^2 - a^2));
sr = theta./rho;

[ng,nc] = size(udg);
nch = 1;
sr_udg = zeros(ng,nch,nc); 

