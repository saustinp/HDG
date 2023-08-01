function eta = get_eta(E, N)
    EN = E/N/1e-21;  % Convert reduced electric field to Td
    log_EN = log10(EN);

    % % Curve fits are only valid for [0.1, 1000] Td
    % if any(log_EN<-1)
    %     error('Input E field out of bounds: lower bound')
    % elseif any(log_EN>3)
    %     error('Input E field out of bounds: upper bound')
    % end
    
    step2 = 0.5*(tanh(1000*(EN-1000))+1);       % This is in linear units
    step_eta = 0.5*(tanh(1000*(log_EN-1.1))+1);
    eta_reg1 = -40.12938813577186 -1.044327759364851*log_EN;
    eta_reg2 = -125.18044468718486 + 168.33587943462615*log_EN -103.64966362614207*log_EN.^2 + 28.456906756759984*log_EN.^3 -2.9427386149316104*log_EN.^4;
    eta_combined_log = (1-step_eta).*eta_reg1 + step_eta.*eta_reg2;
    eta = (1-step2).*power(10, eta_combined_log) + step2.*1.25e-23;
    eta = eta*N;
end