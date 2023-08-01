function alpha = get_alpha(E, N)
    EN = E/N/1e-21;  % Convert reduced electric field to Td
    log_EN = log10(EN);

    % % Curve fits are only valid for [0.1, 1000] Td
    % if any(log_EN<-1)
    %     error('Input E field out of bounds: lower bound')
    % elseif any(log_EN>3)
    %     error('Input E field out of bounds: upper bound')
    % end

    % Ionization/attachment coefficients
    step2 = 0.5*(tanh(1000*(EN-1000))+1);       % This is in linear units
    step_alpha = 0.5*(tanh(1000*(EN-20.09))+1);
    alpha_log = -212.3274758440036 + 304.25126942475583*log_EN -186.28749042736393*log_EN.^2 + 51.50075493521921*log_EN.^3 -5.3618795671332915*log_EN.^4;
    alpha = (1-step2).* step_alpha.*power(10, alpha_log) + step2*1.522e-20;
    alpha = alpha*N;
end