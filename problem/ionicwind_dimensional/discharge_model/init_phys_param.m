function param = init_phys_param()
    % Returns all of the simulation physical parameters that don't depend on the electric field

    % Physical constants
    e = 1.60217663e-19;      % Elementary charge on electron, in C
    epsilon0 = 8.854e-12;    % absolute permittivity of air [C^2/(N*m^2)]
    k_b = 1.380649e-23;      % Boltzmann constant, m2 kg s-2 K-1

    % Static swarm parameters
    mu_p = 2e-4;          % Pos ion mobility [m^2/(Vs)]
    Dp = mu_p*k_b*300/e;     % Einstein relation for ion diffusion coefficient, at 300K. Works out to be 5.1704e-06
    Kep = 2e-13;             % Recombination coeff - pos and neg ions [m^3/s]

    % Same values used for the negative ions, roughly charge independent
    mu_n = mu_p;           % Neg mobility [m^2/(Vs)]
    Dn = Dp;
    Knp = Kep;             % Recombination coeff - pos ions and electrons [m^3/s]

    % Initial condition parameters
    Nmax = 1e16;             % Max number density for initial charge distribution [particles/m^3]
    r0 = 0.0;                % r-pos of emitter tip in reference frame [m]
    z0 = 0.0;                % z-pos of emitter tip in reference frame [m]
    s0 = 40e-6;              % Std deviation of initial charge distribution [m] -> 25e-6 in paper

    % Physical constants and miscellaneous
    Ua = -10e3;              % Emitter potential relative to ground [V]
    gamma = 0.001;           % Secondary electron emission coefficient [1/m]
    E_bd = 3e6;              % Breakdown E field in air [V/m]
    r_tip = 220e-6;          % Tip radius of curvature [m]
    n_ref = epsilon0*E_bd/(e*r_tip);  % Yes, this is redundant and could be recomputed from the above variables. But it saves having to recompute it each time in the functions.
    % ^=7.5357e+17

    P = 101325; % Ambient pressure, Pa
    V = 1; % m^3
    T = 273.15; % K
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