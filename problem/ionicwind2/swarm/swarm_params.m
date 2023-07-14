function params = swarm_params(E, N)
    EN = E/N/1e-21;  % Convert reduced electric field to Td
    log_EN = log10(EN);

    % Curve fits are only valid for [0.1, 1000] Td
    if any(log_EN<-1)
        error('Input E field out of bounds: lower bound')
    elseif any(log_EN>3)
        error('Input E field out of bounds: upper bound')
    end

    % Electron mobility
    mu_e = 10.^24.76713228 * 10.^(-0.34808749.*log_EN);
    
    % % Electron diffusion
    step = 0.5*(tanh(1000*(log_EN-1.9))+1);
    diffusion = (1-step).*(3.477702502838188e23*log_EN + 1.4137092744995724e24)  + step.*(4.6338371114205166e24*log_EN + -6.700654851204797e24);
    
    % Ionization/attachment coefficients
    step_alpha = 0.5*(tanh(1000*(EN-20.09))+1);
    alpha_log = -212.3274758440036 + 304.25126942475583*log_EN -186.28749042736393*log_EN.^2 + 51.50075493521921*log_EN.^3 -5.3618795671332915*log_EN.^4;
    alpha_combined = step_alpha.*power(10, alpha_log);
    
    step_eta = 0.5*(tanh(1000*(log_EN-1.1))+1);
    eta_reg1 = -40.12938813577186 -1.044327759364851*log_EN;
    eta_reg2 = -125.18044468718486 + 168.33587943462615*log_EN -103.64966362614207*log_EN.^2 + 28.456906756759984*log_EN.^3 -2.9427386149316104*log_EN.^4;
    eta_combined_log = (1-step_eta).*eta_reg1 + step_eta.*eta_reg2;
    eta_combined = power(10, eta_combined_log);

    params = {mu_e/N, diffusion/N, alpha_combined*N, eta_combined*N};     % Adjusting by N because coeffs are provided with different units

end
