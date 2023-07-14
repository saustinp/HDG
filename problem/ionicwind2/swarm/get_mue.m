function mu_e = get_mue(E, N)
    EN = E/N/1e-21;  % Convert reduced electric field to Td
    log_EN = log10(EN);

    % % Curve fits are only valid for [0.1, 1000] Td
    % if any(log_EN<-1)
    %     error('Input E field out of bounds: lower bound')
    % elseif any(log_EN>3)
    %     error('Input E field out of bounds: upper bound')
    % end

    % Electron mobility
    mu_e = 10.^24.76713228 * 10.^(-0.34808749.*log_EN);
    mu_e = mu_e/N;
end