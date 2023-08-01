% Script to validate the swarm parameters

% Load in data from BOLSIG swarm paramter file
load bolsig_data.mat BOLSIG_data        % Load swarm parameters from BOLSIG results
E_bolsig = BOLSIG_data(:,1);
mu_bolsig = BOLSIG_data(:,2);
D_bolsig = BOLSIG_data(:,3);
alpha_bolsig = BOLSIG_data(:,4);
eta_bolsig = BOLSIG_data(:,5);

% Data that the curve fits will be loaded into
N = 2.686780111798444e+25;
npts = 200;
E_validate = logspace(-7,4,npts)*N*1e-21;     % This is an electric field vector in [V/m], manufactured from a sweep of the reduced electric field and then multiplied by N and 1e-21

% param = swarm_params(E_validate, N); -> deprecated, old version
% [mue, D, alpha, eta] = param{:};
mue = get_mue(E_validate, N);
D = get_diffusion_e(E_validate, N);
alpha = get_alpha(E_validate, N);
eta = get_eta(E_validate, N);

EN = E_validate/N/1e-21;
figure();
loglog(EN, mue*N, LineWidth=3);  hold on;
loglog(BOLSIG_data(:,1), BOLSIG_data(:,2), LineWidth=3);
title('Electron mobility*N');
xlabel('E/N [Td]');
ylabel('Electron mobility*N');
legend('Curve fit', 'BOLSIG');

figure();
semilogx(EN, D*N, LineWidth=3);  hold on;
loglog(BOLSIG_data(:,1), BOLSIG_data(:,3), LineWidth=3);
title('Electron diffusion coefficient*N');
xlabel('E/N [Td]');
ylabel('Electron diffusion coefficient*N');
legend('Curve fit', 'BOLSIG');

figure();
loglog(EN, alpha/N, LineWidth=3);  hold on;
loglog(EN, eta/N, LineWidth=3);  hold on;
loglog(BOLSIG_data(:,1), BOLSIG_data(:,4), LineWidth=3);  hold on;
loglog(BOLSIG_data(:,1), BOLSIG_data(:,5), LineWidth=3);  hold on;
title('alpha/N, eta/N');
xlabel('E/N [Td]');
ylabel('coefficient/N');
legend('Curve fit: alpha', 'Curve fit: alpha', 'BOLSIG: alpha', 'BOLSIG: alpha');

% figure();
% loglog(EN, alpha/N-eta/N, LineWidth=3);  hold on;
% title('alpha/N - eta/N');
% xlabel('E/N [Td]');
% ylabel('Reduced net production coefficient');
% legend('Reduced net production coefficient');