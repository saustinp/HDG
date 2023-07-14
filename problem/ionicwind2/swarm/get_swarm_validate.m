% Script to validate the swarm parameters

% Load in data from BOLSIG swarm paramter file
load bolsig_data.mat BOLSIG_data        % Load swarm parameters from BOLSIG results
E_bolsig = BOLSIG_data(:,1);
mu_bolsig = BOLSIG_data(:,2);
D_bolsig = BOLSIG_data(:,3);
alpha_bolsig = BOLSIG_data(:,4);
eta_bolsig = BOLSIG_data(:,5);

% Data that the curve fits will be loaded into
npts = 200;
E_validate = logspace(-1,3,npts)*N*1e-21;     % This is an electric field vector in [V/m], manufactured from a sweep of the reduced electric field and then multiplied by N and 1e-21

N = 2.686780111798444e+25;

% E_validate = logspace()
% alpha_vec = zeros(npts);
% eta_vec = zeros(npts);
% D_vec = zeros(npts);
% mue_vec = zeros(npts);
% mun_vec = zeros(npts);

param = get_swarm_params(E, N);
[alpha, eta, beta, D, mue, mup, mun] = swarm{:};

% for i=1:size(E,2)
%     swarm = get_swarm_params(E, N);

%     [alpha, eta, beta, D, mue, mup, mun] = swarm{:};
%     disp(mup);
%     alpha_vec(i) = alpha;
%     eta_vec(i) = eta;
%     D_vec(i) = D;
%     mue_vec(i) = mue;
%     mun_vec(i) = mun;
% end

EN = E_validate/N/1e-21;
figure();
loglog(EN, mue);  hold on;
loglog(BOLSIG_data(:,1), BOLSIG_data(:,2));
title('Electron mobility/N')
xlabel('E/N [Td]')

figure();
semilogx(EN, D, 'black');  hold on;
loglog(BOLSIG_data(:,1), BOLSIG_data(:,3));
title('Electron diffusion coefficient/N')
xlabel('E/N [Td]')

figure();
loglog(EN, alpha);  hold on;
loglog(EN, eta);  hold on;
loglog(BOLSIG_data(:,1), BOLSIG_data(:,4));  hold on;
loglog(BOLSIG_data(:,1), BOLSIG_data(:,5));  hold on;
title('alpha/N, eta/N')
xlabel('E/N [Td]')