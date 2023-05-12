function uex = exactsol(p, param)

K = param(1);
e0 = param(3);
Ec = param(4);
a = param(5);
b = param(6);
J0 = param(7);

d  = J0/(2*pi*e0*K);
c  = a*sqrt(Ec^2 - d);

x = p(:,1,:);
y = p(:,2,:);
r = sqrt(x.^2 + y.^2);

f = sqrt(d*b^2 + c^2);
g = sqrt(d*r.^2 + c^2);

phi = -(g - f + c*log((c+f)./(c+g)) + c*log(r/b));
phix = (x.*(d*r.^2 - a^2*d + Ec^2*a^2 + a*(d*r.^2 - a^2*(- Ec^2 + d)).^(1/2)*(Ec^2 - d)^(1/2)))./(r.^2.*((d*r.^2 - a^2*(- Ec^2 + d)).^(1/2) + a*(Ec^2 - d).^(1/2)));
phiy = (y.*(d*r.^2 - a^2*d + Ec^2*a^2 + a*(d*r.^2 - a^2*(- Ec^2 + d)).^(1/2)*(Ec^2 - d)^(1/2)))./(r.^2.*((d*r.^2 - a^2*(- Ec^2 + d)).^(1/2) + a*(Ec^2 - d).^(1/2)));
 
rho = d./sqrt(d*(r.^2-a^2)+ a^2*Ec^2);
rhox = (d^2*x)./(Ec^2*a^2 + d*(- a^2 + r.^2)).^(3/2);
rhoy = (d^2*y)./(Ec^2*a^2 + d*(- a^2 + r.^2)).^(3/2);
 
uex(:,1,:) = phi;
uex(:,2,:) = rho;
uex(:,3,:) = phix;
uex(:,4,:) = rhox;
uex(:,5,:) = phiy;
uex(:,6,:) = rhoy;


% f = sqrt(J0*b^2/d+c^2);
% g = sqrt(J0*r.^2/d + c^2);
% phi = g - f + c*log((c+f)./(c+g)) + c*log(r/b);
% rho = -(J0/d)./sqrt((J0/d)*(1-(a^2)./(r.^2))+ (a^2*Ec^2)./(r.^2));


return;

syms a b c d Ec x y r

r = sqrt(x^2 + y^2);
c  = a*sqrt(Ec^2 - d);
f = sqrt(d*b^2 + c^2);
g = sqrt(d*r^2 + c^2);
phi = -(g - f + c*log((c+f)./(c+g)) + c*log(r/b));
rho = d./sqrt(d*(r^2-a^2)+ a^2*Ec^2);

phix = simplify(diff(phi,'x'));
phixx = simplify(diff(phix,'x'));
phiy = simplify(diff(phi,'y'));
phiyy = simplify(diff(phiy,'y'));
rhox = simplify(diff(rho,'x'));
rhoy = simplify(diff(rho,'y'));

ddphi = simplify(phixx+phiyy);
simplify(ddphi+rho)
simplify(rho*rho-phix*rhox-phiy*rhoy)


% d  = 2*pi*e0*K;
% c  = a*sqrt(Ec^2 - J0/(d));
% r = sqrt(x.^2 + y.^2);
% f = sqrt(J0*b^2/d+c^2);
% g = sqrt(J0*r.^2/d + c^2);
% phi = g - f + c*log((c+f)./(c+g)) + c*log(r/b);
% rho = -(J0./(2*pi*r*K))./sqrt((J0/d)*(1-(a^2)./(r.^2))+ (a^2*Ec^2)./(r.^2));
% 
% phix = simplify(diff(phi,'x'));
% phixx = simplify(diff(phix,'x'));
% phiy = simplify(diff(phi,'y'));
% phiyy = simplify(diff(phiy,'y'));
% ddphi = simplify(phixx+phiyy);


