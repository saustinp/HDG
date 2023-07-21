function q0 = initq_func_z(dgnodes, param)
    r = dgnodes(:,1,:);
    z = dgnodes(:,2,:);

    % Physics parameters
    Nmax = param{7};
    r0 = param{8};
    z0 = param{9};
    s0 = param{10};

    % Gradient is intentionally negative because q = -grad(u)
    q0 = Nmax* exp(-1/(2*s0.^2) * ((r-r0).^2 + (z-z0).^2)) .* ((z-z0)./s0.^2);
end