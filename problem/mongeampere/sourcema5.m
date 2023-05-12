function [sr, sr_udg] = sourcema5(p,udg,param,time)

theta = param(1);
a1 = param(2);
a2 = param(3);
a = param(4);

qx = udg(:,2);
qy = udg(:,3);
r = sqrt(qx.^2 + qy.^2);
rho = 1 + a1*sech(a2*(r.^2 - a^2));
sr = theta./rho;

[ng,nc] = size(udg);
nch = 1;
sr_udg = zeros(ng,nch,nc); 

tm = (2*a1*a2*theta*sinh(a2*(- a^2 + qx.^2 + qy.^2)))./(a1 + cosh(a2*(- a^2 + qx.^2 + qy.^2))).^2;
sr_udg(:,1,2) = qx.*tm;
sr_udg(:,1,3) = qy.*tm;

