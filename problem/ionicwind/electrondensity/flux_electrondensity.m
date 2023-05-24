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
nch = nc/3;     % 3 components to UDG: (U,QX,QY) -> assuming 2D. For example, in 2D if U has 4 equations, nc will be 12
nq = nc-nch;
nd = nq/nch;

% Physics parameters
% r0 = param{1};
% z0 = param{2};
% s0 = param{3};
% Nmax = param{4};
% e = param{5};
% epsilon0 = param{6};
% Ua = param{7};
% gamma = param{8};
E_bd = param{9};
r_tip = param{10};
% n_ref = param{11};
% N = param{12};
mue_ref = param{13};

% These swarm parameters will eventually be replaced with calls to the curve fit functions.
mu_e = .0378;
mu_p = 2.43e-4;          % mu[3] Pos ion mobility [m^2/(Vs)]
mu_n = 2.7e-4;           % mu[4] Neg mobility [m^2/(Vs)]
De = 0.18;               % Electron diffusion coefficient [m^2/s]
Dn = 0.043e-4;           % Neg diffusion coefficient [m^2/s]
Dp = 0.028e-4;           % Pos ion diffusion coefficient [m^2/s]

De_star = De/(mue_ref*E_bd*r_tip);
Dn_star = Dn/(mue_ref*E_bd*r_tip);
Dp_star = Dp/(mue_ref*E_bd*r_tip);

% Note the u vector is: (u1, u2, u3, u4, q1x, q2x, q3x, q4x, q1y, q2y, q3y, q4y)
%                        1   2   3    4     5       6       7     8     9       10      11    12
% Note the u vector is: (ne, np, nn, phi, dne_dr, dnn_dr, dnp_dr, Er, dne_dz, dnn_dz, dnp_dz, Ez)

% Without the electrostatic equation:
% Note the u vector is: (u1, u2, u3, q1x, q2x, q3x, q1y, q2y, q3y)
%                        1   2   3     4       5       6       7       8       9
% Note the u vector is: (ne, np, nn, dne_dr, dnn_dr, dnp_dr, dne_dz, dnn_dz, dnp_dz)

r = p(:,1);
Ex = p(:,3);    % Ex is -grad(phi)
Ey = p(:,4);

% Read in values from the u vector
ne = udg(:,1);
nn = udg(:,2);
np = udg(:,3);
dne_dr = udg(:,4); % q is -grad(ne)
dnn_dr = udg(:,5); % q is -grad(ne)
dnp_dr = udg(:,6);
dne_dz = udg(:,7);
dnn_dz = udg(:,8);
dnp_dz = udg(:,9);

% Compute convective velocities
cr_e = -(mu_e/mue_ref)*Ex;
cz_e = -(mu_e/mue_ref)*Ey;
cr_n = -(mu_n/mue_ref)*Ex;
cz_n = -(mu_n/mue_ref)*Ey;
cr_p = (mu_p/mue_ref)*Ex;
cz_p = (mu_p/mue_ref)*Ey;

% For the 3-species problem (no poisson yet), nch=3, nd=2
f = zeros(ng,nch,nd);
f(:,1,1) = cr_e.*(r.*ne) + De_star*(r.*dne_dr);    % fr, q is -grad(ne)
f(:,1,2) = cz_e.*(r.*ne) + De_star*(r.*dne_dz);    % fz, q is -grad(ne)
f(:,2,1) = cr_n.*(r.*nn) + Dn_star*(r.*dnn_dr);
f(:,2,2) = cz_n.*(r.*nn) + Dn_star*(r.*dnn_dz);
f(:,3,1) = cr_p.*(r.*np) + Dp_star*(r.*dnp_dr);
f(:,3,2) = cz_p.*(r.*np) + Dp_star*(r.*dnp_dz);

% Jacobian of flux with respect to UDG. Thesea are the only nonzero components, the rest (cross terms etc) are zero. For example dfne_r_d(dne_dz) - f_udg(:,1,1,3)
% This will be more complicated when the electrostatic equation is added
f_udg = zeros(ng,nch,nd,nc);
f_udg(:,1,1,1) = cr_e.*r;       % df(ne)_r_d(ne)
f_udg(:,1,1,4) = De_star*r;     % df(ne)_r_d(dne_dr)
f_udg(:,1,2,1) = cz_e.*r;       % df(ne)_z_d(ne)
f_udg(:,1,2,7) = De_star*r;     % df(ne)_z_d(dne_dz)

f_udg(:,2,1,2) = cr_n.*r;       % df(nn)_r_d(nn)
f_udg(:,2,1,5) = Dn_star*r;     % df(nn)_r_d(dnn_dr)
f_udg(:,2,2,2) = cz_n.*r;       % df(nn)_z_d(nn)
f_udg(:,2,2,8) = Dn_star*r;     % df(nn)_z_d(dnn_dz)

f_udg(:,3,1,3) = cr_p.*r;       % df(np)_r_d(np)
f_udg(:,3,1,6) = Dp_star*r;     % df(np)_r_d(dnp_dr)
f_udg(:,3,2,3) = cz_p.*r;       % df(np)_z_d(np)
f_udg(:,3,2,9) = Dp_star*r;     % df(np)_z_d(dnp_dz)