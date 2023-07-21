function diffusion_coefficent = get_diffusion_e(E, N)
    EN = E/N/1e-21;  % Convert reduced electric field to Td
    log_EN = log10(EN);

    % % Curve fits are only valid for [0.1, 1000] Td
    % if any(log_EN<-1)
    %     error('Input E field out of bounds: lower bound')
    % elseif any(log_EN>3)
    %     error('Input E field out of bounds: upper bound')
    % end
    
    % Electron diffusion
    step = 0.5*(tanh(1000*(log_EN-1.9))+1);
    diffusion_coefficent = (1-step).*(3.477702502838188e23*log_EN + 1.4137092744995724e24)  + step.*(4.6338371114205166e24*log_EN + -6.700654851204797e24);
    diffusion_coefficent = diffusion_coefficent/N;
end