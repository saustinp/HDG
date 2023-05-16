function u0 = initu_func(dgnodes, physparam)
    % size(dgnodes)
    r = dgnodes(:,1,:);
    z = dgnodes(:,2,:);

    % Physics parameters
    r0 = physparam{1};
    z0 = physparam{2};
    s0 = physparam{3};
    Nmax = physparam{4};
    e = physparam{5};
    epsilon0 = physparam{6};
    Ua = physparam{7};
    gamma = physparam{8};
    E_bd = physparam{9};
    r_tip = physparam{10};
    n_ref = physparam{11};
    N = physparam{12};
    mue_ref = physparam{13};
    D_star = physparam{14};

    % Nondimensionalization
    r0_star = r0/r_tip;
    z0_star = z0/r_tip;
    sig0_star = s0/r_tip;
    n_max_star = Nmax/n_ref;

    % "star" indicating it is a nondimensional quantity
    n0_star = n_max_star* exp(-1/(2*sig0_star.^2) * ((r-r0_star).^2 + (z-z0_star).^2));

    u0 = n0_star;    % Only set the seed particles for the electrons and positives
end