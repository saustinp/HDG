function [f,f_udg] = flux_electrondensity(p,udg,param,time)
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

% Physics parameters
% r0 = param{1};
% z0 = param{2};
% s0 = param{3};
% Nmax = param{4};
% e = param{5};
% epsilon0 = param{6};
% Ua = param{7};
% gamma = param{8};
% E_bd = param{9};
% r_tip = param{10};
% n_ref = param{11};
% N = param{12};
mue_ref = param{13};
D_star = param{14};
mu_e = .0378;

r = p(:,1);
Ex = p(:,3);    % Ex is -grad(phi)
Ey = p(:,4);
ne = udg(:,1);
dne_dr = udg(:,2); % q is -grad(ne)
dne_dz = udg(:,3); % q is -grad(ne)
c1 = -(mu_e/mue_ref)*Ex;
c2 = -(mu_e/mue_ref)*Ey;

% f = r*(-(mu_e/mue_ref)*[Ex, Ey].*ne + D_star*[dne_dr dne_dz]);    From Exasim

f = zeros(ng,nch,nd);
f(:,1,1) = c1.*(r.*ne) + D_star*(r.*dne_dr);    % fr, q is -grad(ne)
f(:,1,2) = c2.*(r.*ne) + D_star*(r.*dne_dz);    % fz, q is -grad(ne)

f_udg = zeros(ng,nch,nd,nc);
f_udg(:,1,1,1) = c1.*r;    
f_udg(:,1,1,2) = D_star*r;
f_udg(:,1,2,1) = c2.*r;
f_udg(:,1,2,3) = D_star*r;