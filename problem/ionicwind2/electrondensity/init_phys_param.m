function param = init_phys_param()
    % Returns all of the simulation physical parameters that don't depend on the electric field

    % Static swarm parameters
    mu_p = 2e-4;          % Pos ion mobility [m^2/(Vs)]
    mu_n = 2e-4;           % Neg mobility [m^2/(Vs)]
    Dn = 5.1704e-6;           % Neg diffusion coefficient [m^2/s]
    Dp = 5.1704e-6;           % Pos ion diffusion coefficient [m^2/s]
    Kep = 2e-13;             % Recombination coeff - pos and neg ions [m^3/s]
    Knp = 2e-13;             % Recombination coeff - pos ions and electrons [m^3/s]

    % Initial condition parameters
    Nmax = 1e16;             % Max number density for initial charge distribution [particles/m^3]
    r0 = 0.0;                % r-pos of emitter tip in reference frame [m]
    z0 = 0.0;                % z-pos of emitter tip in reference frame [m]
    s0 = 25e-6;              % Std deviation of initial charge distribution [m] -> 25e-6 in paper
    e = 1.6022e-19;          % Charge on electron [C]

    % Physical constants and miscellaneous
    epsilon0 = 8.854e-12;    % absolute permittivity of air [C^2/(N*m^2)]
    Ua = -10e3;              % Emitter potential relative to ground [V]
    gamma = 0.001;           % Secondary electron emission coefficient [1/m]
    E_bd = 3e6;              % Breakdown E field in air [V/m]
    r_tip = 220e-6;          % Tip radius of curvature [m]
    n_ref = epsilon0*E_bd/(e*r_tip);  % Yes, this is redundant and could be recomputed from the above variables. But it saves having to recompute it each time in the functions.
    % ^=7.5357e+17

    P = 101325; % Ambient pressure, Pa
    V = 1; % m^3
    T = 273.15; % K
    k_b = 1.380649e-23; % m2 kg s-2 K-1
    N = P*V/(k_b*T);    % Neutral number density
    mue_ref = 0.04266918357567234;   % m^2/(V*s)    Taken at 3e6 V/m

    %          1     2    3  4    5    6    7     8  9   10  11   12      13   14     15    16      17  18   19
    param = {mu_p, mu_n, Dn, Dp, Kep, Knp, Nmax, r0, z0, s0, e, epsilon0, Ua, gamma, E_bd, r_tip, n_ref, N, mue_ref};

    % Deprecated: output from swarm params script:
    % N: 2.4614924955148245e+25 particles/m^3
    % At 130 Td:
    % alpha static: 1036.2883406117412
    % eta static: 975.9817744716279
    % diffusion static: 0.1185053379327856
    % mobility static: 0.04627273908311374
end

% Reference list of physics parameters
% mu_p = param(1);
% mu_n = param(2);
% Dn = param(3);
% Dp = param(4);
% Kep = param(5);
% Knp = param(6);
% Nmax = param(7);
% r0 = param(8);
% z0 = param(9);
% s0 = param(10);
% e = param(11);
% epsilon0 = param(12);
% Ua = param(13);
% gamma = param(14);
% E_bd = param(15);
% r_tip = param(16);
% n_ref = param(17);
% N = param(18);
% mue_ref = param(19);
% tau = param(end);