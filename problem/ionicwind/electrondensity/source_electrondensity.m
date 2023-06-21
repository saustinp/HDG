function [sr,sr_udg] = source_electrondensity(p,udg,param,time)

[ng,nc] = size(udg);
nch = nc/3;     % 3 components to UDG: (U,QX,QY) -> assuming 2D. For example, in 2D if U has 4 equations, nc will be 12

% Physics parameters
% r0 = param{1};
% z0 = param{2};
% s0 = param{3};
% Nmax = param{4};
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

% Output from swarm param script:
% N: 2.4614924955148245e+25 particles/m^3
% At 130 Td:
% alpha static: 1036.2883406117412
% eta static: 975.9817744716279
% diffusion static: 0.1185053379327856
% mobility static: 0.04627273908311374


% alpha = 1036.28;
% eta = 975.98;

% beta = 2e-13;
% mue = .0378;
% mup = 2.34e-4;
% mun = -2.7e-4;

r = p(:,1);

sr = zeros(ng,nch);
sr_udg = zeros(ng,nch,nc);

Er_0 = p(:,3);    % Er is -grad(phi)
Ez_0 = p(:,4);

% Read in values from the u vector
ne = udg(:,1);
nn = udg(:,2);
np = udg(:,3);

% step_ne = 0.5*(tanh(1000*(ne))+1);
% step_np = 0.5*(tanh(1000*(np))+1);
% step_nn = 0.5*(tanh(1000*(nn))+1);
% ne = (1-step_ne)*1e-5 + step*ne
% np = (1-step_np)*1e-5 + step*np
% nn = (1-step_nn)*1e-5 + step*nn

f = @(x, alpha) x.*(atan(alpha*x)/pi + 1/2) + 1/2 - atan(alpha)/pi;
% f_prime = @(x) atan(1000.*x)/pi + (1000.*x)/(pi*(1e6.*x.^2 + 1)) + 1/2;

ne = f(ne, 1000);
np = f(np, 1000);
nn = f(nn, 1000);

phi = udg(:,4);
dne_dr = udg(:,5); % q is -grad(ne)
dnn_dr = udg(:,6); % q is -grad(ne)
dnp_dr = udg(:,7);
Er = udg(:,8) + Er_0;
dne_dz = udg(:,9);
dnn_dz = udg(:,10);
dnp_dz = udg(:,11);
Ez = udg(:,12) + Ez_0;

normE = sqrt(Er.^2+Ez.^2);

swarm = get_swarm_params(normE, N);
alpha = swarm(:,1);
eta = swarm(:,2);
beta = swarm(:,3);
De = swarm(:,4);
mue = swarm(:,5);
mup = swarm(:,6);
mun = swarm(:,7);

k_ep = beta;
k_np = beta;


% Phase 1: no source terms - this is the end of the function


% % Phase 2: source terms for only the first equation
% % se
% sr(:,1) = (alpha-eta)*(mue/mue_ref)*r_tip*ne.*normE;

% % Jacobian of sr with respect to u -> (U and Q). Only dse_dne term is non-zero
% % dse_dne
% sr_udg(:,1,1) = (mue.*r_tip.*(alpha - eta).*(Er.^2 + Ez.^2).^(1./2))./mue_ref;

% Note that a factor of 'r' was added for axisymmetric
% Phase 3: Source terms for all the equations
sr(:,1) = r.*((alpha-eta).*(mue/mue_ref).*r_tip.*ne.*normE - k_ep.*epsilon0/(e.*mue_ref).*ne.*np);                                  % ne
sr(:,2) = r.*(eta.*(mue/mue_ref).*r_tip.*ne.*normE - k_np.*epsilon0/(e.*mue_ref).*nn.*np);                                          % nn
sr(:,3) = r.*(alpha.*(mue/mue_ref).*r_tip.*ne.*normE - np.*epsilon0./(e.*mue_ref).*(k_np.*nn + k_ep.*ne));                            % np
sr(:,4) = -r.*(ne + nn - np);                                                                                                 % phi

% Jacobian of sr with respect to u -> (U and Q)
sr_udg(:,1,1) = r.*((mue.*r_tip.*(alpha - eta).*(Er.^2 + Ez.^2).^(1./2))./mue_ref - (epsilon0.*k_ep.*np)./(e.*mue_ref));     % dse_dne
sr_udg(:,1,2) = 0;                                                                                                          % dse_dnn
sr_udg(:,1,3) = r.*(-(epsilon0.*k_ep.*ne)./(e.*mue_ref));                                                                    % dse_dnp
sr_udg(:,1,8) = r.*((Er.*mue.*ne.*r_tip.*(alpha - eta))./(mue_ref.*(Er.^2 + Ez.^2).^(1./2)));                                % dse_dEr
sr_udg(:,1,12) = r.*((Ez.*mue.*ne.*r_tip.*(alpha - eta))./(mue_ref.*(Er.^2 + Ez.^2).^(1./2)));                               % dse_dEz

sr_udg(:,2,1) = r.*((eta.*mue.*r_tip.*(Er.^2 + Ez.^2).^(1./2))./mue_ref);                                                    % dsn_dne
sr_udg(:,2,2) = r.*(-(epsilon0.*k_np.*np)./(e.*mue_ref));                                                                    % dsn_dnn
sr_udg(:,2,3) = r.*(-(epsilon0.*k_np.*nn)./(e.*mue_ref));                                                                    % dsn_dnp
sr_udg(:,2,8) =  r.*((Er.*eta.*mue.*ne.*r_tip)./(mue_ref.*(Er.^2 + Ez.^2).^(1./2)));                                          % dsn_dEr
sr_udg(:,2,12) = r.*((Ez.*eta.*mue.*ne.*r_tip)./(mue_ref.*(Er.^2 + Ez.^2).^(1./2)));                                         % dsn_dEz

sr_udg(:,3,1) = r.*((alpha.*mue.*r_tip.*(Er.^2 + Ez.^2).^(1./2))./mue_ref - (epsilon0.*k_ep.*np)./(e.*mue_ref));             % dsp_dne
sr_udg(:,3,2) = r.*(-(epsilon0.*k_np.*np)./(e.*mue_ref));                                                                    % dsp_dnn
sr_udg(:,3,3) = r.*(-(epsilon0.*(k_ep.*ne + k_np.*nn))./(e.*mue_ref));                                                       % dsp_dnp
sr_udg(:,3,8) =  r.*((Er.*alpha.*mue.*ne.*r_tip)./(mue_ref.*(Er.^2 + Ez.^2).^(1./2)));                                        % dsp_dEr
sr_udg(:,3,12) = r.*((Ez.*alpha.*mue.*ne.*r_tip)./(mue_ref.*(Er.^2 + Ez.^2).^(1./2)));                                       % dsp_dEz

sr_udg(:,4,1) = -r;                                                                                                      % dsphi_dne
sr_udg(:,4,2) = -r;                                                                                                      % dsphi_dnn
sr_udg(:,4,3) = r;                                                                                                     % dsphi_dnp
