% N = 2.686780111798444e+25;
% load bolsig_data.mat BOLSIG_data        % Load swarm parameters from BOLSIG results      

% swarm = get_swarm_params(130*2.6868e+04, N);
% disp(swarm)

% 
% alpha_vec = zeros(100);
% eta_vec = zeros(100);
% D_vec = zeros(100);
% mue_vec = zeros(100);
% mun_vec = zeros(100);
% 
% 
% for i=1:size(E,2)
%     swarm = get_swarm_params(E(i), N);
%     swarm = num2cell(swarm);
%     [alpha, eta, beta, D, mue, mup, mun] = swarm{:};
%     disp(mup);
%     alpha_vec(i) = alpha;
%     eta_vec(i) = eta;
%     D_vec(i) = D;
%     mue_vec(i) = mue;
%     mun_vec(i) = mun;
% end
% 
% figure;
% plot(E, alpha_vec);
% title('alpha');
% figure;
% plot(E, eta_vec);
% figure;
% plot(E, D_vec);
% figure;
% plot(E, mue_vec);
% figure;
% plot(E, mun_vec);
% disp(beta)



% EN = logspace(-1,3,200);
% 
% log_EN = log10(EN);
% 
% % Electron mobility
% mu_e = 10.^24.76713228 * 10.^(-0.34808749.*log_EN);
% figure();
% loglog(EN, mu_e);  hold on;
% loglog(BOLSIG_data(:,1), BOLSIG_data(:,2));
% title('Electron mobility/N')
% xlabel('E/N [Td]')
% 
% % % Electron diffusion
% step = 0.5*(tanh(1000*(log_EN-1.9))+1);
% diffusion = (1-step).*(3.477702502838188e23*log_EN + 1.4137092744995724e24)  + step.*(4.6338371114205166e24*log_EN + -6.700654851204797e24);
% figure();
% semilogx(EN, diffusion_combined, 'black');  hold on;
% loglog(BOLSIG_data(:,1), BOLSIG_data(:,3));
% title('Electron diffusion coefficient/N')
% xlabel('E/N [Td]')
% 
% % Ionization/attachment coefficients
% step_alpha = 0.5*(tanh(1000*(EN-20.09))+1);
% alpha_log = -212.3274758440036 + 304.25126942475583*log_EN -186.28749042736393*log_EN.^2 + 51.50075493521921*log_EN.^3 -5.3618795671332915*log_EN.^4;
% alpha_combined = step_alpha.*power(10, alpha_log);
% 
% step_eta = 0.5*(tanh(1000*(log_EN-1.1))+1);
% eta_reg1 = -40.12938813577186 -1.044327759364851*log_EN;
% eta_reg2 = -125.18044468718486 + 168.33587943462615*log_EN -103.64966362614207*log_EN.^2 + 28.456906756759984*log_EN.^3 -2.9427386149316104*log_EN.^4;
% eta_combined_log = (1-step_eta).*eta_reg1 + step_eta.*eta_reg2;
% eta_combined = power(10, eta_combined_log);
% 
% figure();
% loglog(EN, alpha_combined);  hold on;
% loglog(EN, eta_combined);  hold on;
% loglog(BOLSIG_data(:,1), BOLSIG_data(:,4));  hold on;
% loglog(BOLSIG_data(:,1), BOLSIG_data(:,5));  hold on;
% title('alpha/N, eta/N')
% xlabel('E/N [Td]')


function params = get_swarm_params(E, N)

    EN = E/N;
    EN = EN/1e-21;  % Convert to Td
    cm2m = 0.01;

    log_EN = log10(EN);

    % Electron mobility
    mu_e = 10.^24.76713228 * 10.^(-0.34808749.*log_EN);
    % figure();
    % loglog(EN, mu_e);  hold on;
    % loglog(BOLSIG_data(:,1), BOLSIG_data(:,2));
    % title('Electron mobility/N')
    % xlabel('E/N [Td]')
    
    % % Electron diffusion
    step = 0.5*(tanh(1000*(log_EN-1.9))+1);
    diffusion = (1-step).*(3.477702502838188e23*log_EN + 1.4137092744995724e24)  + step.*(4.6338371114205166e24*log_EN + -6.700654851204797e24);
    % figure();
    % semilogx(EN, diffusion_combined, 'black');  hold on;
    % loglog(BOLSIG_data(:,1), BOLSIG_data(:,3));
    % title('Electron diffusion coefficient/N')
    % xlabel('E/N [Td]')
    
    % Ionization/attachment coefficients
    step_alpha = 0.5*(tanh(1000*(EN-20.09))+1);
    alpha_log = -212.3274758440036 + 304.25126942475583*log_EN -186.28749042736393*log_EN.^2 + 51.50075493521921*log_EN.^3 -5.3618795671332915*log_EN.^4;
    alpha_combined = step_alpha.*power(10, alpha_log);
    
    step_eta = 0.5*(tanh(1000*(log_EN-1.1))+1);
    eta_reg1 = -40.12938813577186 -1.044327759364851*log_EN;
    eta_reg2 = -125.18044468718486 + 168.33587943462615*log_EN -103.64966362614207*log_EN.^2 + 28.456906756759984*log_EN.^3 -2.9427386149316104*log_EN.^4;
    eta_combined_log = (1-step_eta).*eta_reg1 + step_eta.*eta_reg2;
    eta_combined = power(10, eta_combined_log);
    
    % figure();
    % loglog(EN, alpha_combined);  hold on;
    % loglog(EN, eta_combined);  hold on;
    % loglog(BOLSIG_data(:,1), BOLSIG_data(:,4));  hold on;
    % loglog(BOLSIG_data(:,1), BOLSIG_data(:,5));  hold on;
    % title('alpha/N, eta/N')
    % xlabel('E/N [Td]')

    % Recombination coefficient
    beta = ones(size(EN))*2e-7*cm2m^3;

    % Positive ion mobility
    % Assuming that P0/P is approximately 1 and that the units for the E field are V/cm
    mun = ones(size(EN))*2.34*cm2m^2;
    mup = ones(size(EN))*2.34*cm2m^2;

    params = [alpha_combined*N, eta_combined*N, beta, diffusion/N, mu_e/N, mup, mun];     % Adjusting by N because coeffs are provided with different units
end