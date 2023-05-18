function [sr,sr_udg] = source_electrondensity(p,udg,param,time)

[ng,nc] = size(udg);
nch = 1;

% Physics parameters
r0 = param{1};
z0 = param{2};
s0 = param{3};
Nmax = param{4};
e = param{5};
epsilon0 = param{6};
Ua = param{7};
gamma = param{8};
E_bd = param{9};
r_tip = param{10};
n_ref = param{11};
N = param{12};
mue_ref = param{13};
D_star = param{14};
Kep = param{15};
Knp = param{16};
mu_p = param{17};
mu_n = param{18};

alpha = 1.19e-21;
eta = 2.28e-19;
D = .16;    % Assuming it's supposed to be positive for now, will need to verify the swarm param plots later
beta = 2e-13;
mue = .0378;
mup = 2.34e-4;
mun = -2.7e-4;

Ex = p(:,3);    % Ex is -grad(phi)
Ey = p(:,4);

% normE = sqrt(Ex.^2 + Ey.^2);
% sr = (alpha-eta)*(mue/mue_ref)*r_tip*ne.*normE - Kep*epsilon0/(e*mue_ref)*ne*np;
% sr_udg = ;

sr = zeros(ng,nch);
sr_udg = zeros(ng,nch,nc);