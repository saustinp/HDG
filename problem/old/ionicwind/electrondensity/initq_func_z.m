function q0 = initq_func_z(dgnodes)
    r = dgnodes(:,1,:);
    z = dgnodes(:,2,:);

    physparam = get_phys_param();

    % Physics parameters
    Nmax = physparam{7};
    r0 = physparam{8};
    z0 = physparam{9};
    s0 = physparam{10};
    r_tip = physparam{16};

    % Nondimensionalization
    r0_star = r0/r_tip;
    z0_star = z0/r_tip;
    sig0_star = s0/r_tip;
    n_max_star = Nmax/n_ref;

    % Gradient is intentionally negative because q = -grad(u)
    n0_star = n_max_star* exp(-1/(2*sig0_star.^2) * ((r-r0_star).^2 + (z-z0_star).^2)) .* (z-z0_star)./sig0_star.^2;

    q0 = n0_star;    % Only set the seed particles for the electrons and positives
end