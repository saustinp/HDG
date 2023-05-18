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
k_ep = param{15};
k_np = param{16};
mu_p = param{17};
mu_n = param{18};

alpha = 1.19e-21;
eta = 2.28e-19;
D = 0.11736375131509072;    % Assuming it's supposed to be positive for now, will need to verify the swarm param plots later
beta = 2e-13;
mue = .0378;
mup = 2.34e-4;
mun = -2.7e-4;

ne = udg(:,1);
Er = p(:,3);    % Ex is -grad(phi)
Ez = p(:,4);
normE = sqrt(Er.^2+Ez.^2);

sr = zeros(ng,nch);
sr_udg = zeros(ng,nch,nc);

% Phase 1: no source terms - this is the end of the function


% Phase 2: source terms for only the first equation
% se
sr(:,1) = (alpha-eta)*(mue/mue_ref)*r_tip*ne.*normE;

% Jacobian of sr with respect to u -> (U and Q). Only dse_dne term is non-zero
% dse_dne
sr_udg(:,1,1) = (mue.*r_tip.*(alpha - eta).*(Er.^2 + Ez.^2).^(1./2))./mue_ref;




% % Phase 3: Source terms for all the equations
% % ne
% sr(:,1) = (alpha-eta)*(mue/mue_ref)*r_tip*ne.*normE - k_ep*epsilon0/(e*mue_ref)*ne*np;
% % np
% sr(:,2) = alpha*(mue/mue_ref)*r_tip*ne.*normE - np*epsilon0/(e*mue_ref)*(k_np*nn + k_ep*ne);
% % nn
% sr(:,3) = eta*(mue/mue_ref)*r_tip*ne.*normE - k_np*epsilon0/(e*mue_ref)*nn*np;
% % phi
% sr(:,4) = ne + nn - np;

% % Jacobian of sr with respect to u -> (U and Q)
% % dse_dne
% sr_udg(:,1,1) = (mue.*r_tip.*(alpha - eta).*(Er.^2 + Ez.^2).^(1./2))./mue_ref - (epsilon0.*k_ep.*np)./(e.*mue_ref);
% % dse_dnp
% sr_udg(:,1,2) = -(epsilon0.*k_ep.*ne)./(e.*mue_ref);
% % dse_dnn
% sr_udg(:,1,3) = 0;
% % dse_dEr
% sr_udg(:,1,8) = (Er.*mue.*ne.*r_tip.*(alpha - eta))./(mue_ref.*(Er.^2 + Ez.^2).^(1./2));
% % dse_dEz
% sr_udg(:,1,12) = (Ez.*mue.*ne.*r_tip.*(alpha - eta))./(mue_ref.*(Er.^2 + Ez.^2).^(1./2));

% % dsp_dne
% sr_udg(:,2,1) = (alpha.*mue.*r_tip.*(Er.^2 + Ez.^2).^(1./2))./mue_ref - (epsilon0.*k_ep.*np)./(e.*mue_ref);
% % dsp_dnp
% sr_udg(:,2,2) = -(epsilon0.*(k_ep.*ne + k_np.*nn))./(e.*mue_ref);
% % dsp_dnn
% sr_udg(:,2,3) = -(epsilon0.*k_np.*np)./(e.*mue_ref);
% % dsp_dEr
% sr_udg(:,2,8) = (Er.*alpha.*mue.*ne.*r_tip)./(mue_ref.*(Er.^2 + Ez.^2).^(1./2));
% % dsp_dEz
% sr_udg(:,2,12) = (Ez.*alpha.*mue.*ne.*r_tip)./(mue_ref.*(Er.^2 + Ez.^2).^(1./2));

% % dsn_dne
% sr_udg(:,3,1) = (eta.*mue.*r_tip.*(Er.^2 + Ez.^2).^(1./2))./mue_ref;
% % dsn_dnp
% sr_udg(:,3,2) = -(epsilon0.*k_np.*nn)./(e.*mue_ref);
% % dsn_dnn
% sr_udg(:,3,3) = -(epsilon0.*k_np.*np)./(e.*mue_ref);
% % dsn_dEr
% sr_udg(:,3,8) = (Er.*eta.*mue.*ne.*r_tip)./(mue_ref.*(Er.^2 + Ez.^2).^(1./2));
% % dsn_dEz
% sr_udg(:,3,12) = (Ez.*eta.*mue.*ne.*r_tip)./(mue_ref.*(Er.^2 + Ez.^2).^(1./2));

% % dsphi_dne
% sr_udg(:,4,1) = 1;
% % dsphi_dnp
% sr_udg(:,4,2) = -1;
% % dsphi_dnn
% sr_udg(:,4,2) = 1;
