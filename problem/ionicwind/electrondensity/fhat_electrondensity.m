function [fh,fh_udg,fh_uh] = fhat_electrondensity(nl,p,udg,uh,param,time,signe)
%FHAT flux function
%   [fh,fhu,fhq,fhm] = fhat(nl,p,u,q,m,param)
%
%      NL(N,ND)              Normal N points
%      P(N,ND)               Coordinates for N points
%      U(N,NC)               Unknown vector for N points with NC components
%      Q(N,NC,ND)            Flux vector for N points with NC components in the
%                            coordinate directions
%      M(N,NC)               Hybrid unkowns for N points with NC components
%      PARAM                 Parameter list
%      FH(N,NC):              Volume flux at N points
%      FHU(N,NC,NC):         Jacobian of the flux flux vector w.r.t. U
%      FHQ(N,NC,NC,ND):      Jacobian of the flux flux vector w.r.t. Q
%      FHQ(N,NC,NC):         Jacobian of the flux flux vector w.r.t. Q

[ng,nc] = size(udg);
nch = 1;
% nq = nc-nch;
% nd = nq;

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
tau   = param{end};
mu_e = .0378;

r = p(:,1);
Ex = p(:,3);    % Ex is -grad(phi)
Ey = p(:,4);
ne = udg(:,1);
dne_dr = udg(:,2); % q is -grad(ne)
dne_dz = udg(:,3); % q is -grad(ne)
c1 = -(mu_e/mue_ref)*Ex;
c2 = -(mu_e/mue_ref)*Ey;

% ng x nch
fh = r.*((c1.*uh+D_star*dne_dr).*nl(:,1) + (c2.*uh+D_star*dne_dz).*nl(:,2)) + tau.*(ne-uh);        % CHANGE ME

fh_udg = zeros(ng,nch,nc);
fh_udg(:,1,1) = tau;
fh_udg(:,1,2) = D_star*r.*nl(:,1);
fh_udg(:,1,3) = D_star*r.*nl(:,2);

fh_uh = zeros(ng,nch,nch);
fh_uh(:,1,1) = c1.*r.*nl(:,1) + c2.*r.*nl(:,2) - tau;