function [f,f_udg] = flux(p,udg,param,time)
%FLUX Volume flux function
%   [f,fu,fq] = flux(p,u,q,param)
%
%      P(N,ND)              Coordinates for N points
%      U(N,NC)              Unknown vector for N points with NC components
%      Q(N,NC,ND)           Flux vector for N points with NC components in the
%                           coordinate directions
%      PARAM                Parameter list
%      F(N,NC,ND):          Volume flux at N points
%      FU(N,NC,ND,NC):      Jacobian of the flux flux vector w.r.t. U
%      FQ(N,NC,ND,NC,ND):   Jacobian of the flux flux vector w.r.t. Q

[ng,nc] = size(udg);
nch = 1;
nq = nc-nch;
nd = nq;

kappa = param{1};
c     = param{2};
%nu = param{3};
%h = param{4};

f = zeros(ng,nch,nd);
f_u = zeros(ng,nch,nd,nch);
f_q = zeros(ng,nch,nd,nd);

u = udg(:,1);
qx = udg(:,2);
qy = udg(:,3);
s = p(:,nd+1);

f(:,1,1) = c(1)*(u.*u) + (s+kappa).*qx;
f(:,1,2) = u + c(2)*(u.*u) + 0*(s+kappa).*qy;

f_u(:,1,1,1) = 2*c(1)*u;
f_u(:,1,2,1) = 1.0 + 2*c(2)*u;
f_q(:,1,1,1) = s + kappa;
f_q(:,1,2,2) = 0*(s + kappa);

if nd==3
    qz = udg(:,4);
    f(:,1,3) = c(3)*(u.*u) + kappa*qz;
    f_u(:,1,3,1) = 2*c(3)*u;
    f_q(:,1,3,3) = kappa;
end

f_udg = cat(4,f_u,f_q);


