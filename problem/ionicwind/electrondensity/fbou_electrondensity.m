function [fh,fh_udg,fh_uh] = fbou_electrondensity(ib,ui,nl,p,udg,uh,param,time)
%FHAT flux function
%   [fh,fhu,fhq,fhm] = fhat(nl,p,u,q,m,param)
%
%      NL(N,ND)              Normal N points
%      P(N,ND)               Coordinates for N points
%      U(N,NC)               Unknown vector for N points with NC components
%      Q(N,NC,ND)            Flux vector for N points with NC components in the
%                            coordinate directions
%      M(N,NC)               Hybrid unkowns for N points with NC components
%      PARAM                 Parameter list
%      FH(N,NC):              Volume flux at N points
%      FHU(N,NC,NC):         Jacobian of the flux flux vector w.r.t. U
%      FHQ(N,NC,NC,ND):      Jacobian of the flux flux vector w.r.t. Q
%      FHQ(N,NC,NC):         Jacobian of the flux flux vector w.r.t. Q

[ng,nc] = size(udg);
nch = nc/3;     % 3 components to UDG: (U,QX,QY) -> assuming 2D. For example, in 2D if U has 4 equations, nc will be 12
% nq = nc-nch;

% Physics parameters
% r0 = param{1};
% z0 = param{2};
% s0 = param{3};
% Nmax = param{4};
% e = param{5};
% epsilon0 = param{6};
% Ua = param{7};
% gamma = param{8};
E_bd = param{9};
r_tip = param{10};
% n_ref = param{11};
% N = param{12};
mue_ref = param{13};

% These swarm parameters will eventually be replaced with calls to the curve fit functions.
% mu_e = .0378;
De = 0.18;               % Electron diffusion coefficient [m^2/s]
Dn = 0.043e-4;           % Neg diffusion coefficient [m^2/s]
Dp = 0.028e-4;           % Pos ion diffusion coefficient [m^2/s]

De_star = De/(mue_ref*E_bd*r_tip);
Dn_star = Dn/(mue_ref*E_bd*r_tip);
Dp_star = Dp/(mue_ref*E_bd*r_tip);
tau = param{end};

r = p(:,1);
Ex_0 = p(:,3);    % Ex is -grad(phi)
Ey_0 = p(:,4);

% Read in values from the u vector
ne = udg(:,1);
nn = udg(:,2);
np = udg(:,3);
phi = udg(:,4);
dne_dr = udg(:,5); % q is -grad(ne)
dnn_dr = udg(:,6); % q is -grad(ne)
dnp_dr = udg(:,7);
Ex = udg(:,8) + Ex_0;
dne_dz = udg(:,9);
dnn_dz = udg(:,10);
dnp_dz = udg(:,11);
Ey = udg(:,12) + Ey_0;
Ex_prime = udg(:,8);
Ey_prime = udg(:,12);

switch ib
    case 1  % Axisymmetry boundary
        fh = zeros(ng,nch);
        fh_udg = zeros(ng,nch,nc);
        fh_uh = zeros(ng,nch,nch);

        % All equations: symmetry BC
        r = 1.0;
        D_star = 1.0;

        % The following is just copied from the fhat_electrondensity function and then the convective component taken out with all the same D_star

        % ng x nch, same size as uh
        % This is the jacobian of the normal flux f_hat w.r.t UDG. UDG is [u, uq1x, q2x], and does _not_ involve uh.
        fh(:,1) = r.*(D_star*dne_dr.*nl(:,1) + D_star*dne_dz.*nl(:,2)) + tau.*(ne-uh(:,1));
        fh(:,2) = r.*(D_star*dnn_dr.*nl(:,1) + D_star*dnn_dz.*nl(:,2)) + tau.*(nn-uh(:,2));
        fh(:,3) = r.*(D_star*dnp_dr.*nl(:,1) + D_star*dnp_dz.*nl(:,2)) + tau.*(np-uh(:,3));
        fh(:,4) = r.*(Ex_prime.*nl(:,1) + Ey_prime.*nl(:,2)) + tau.*(phi-uh(:,4));

        % Jacobian of f_hat with respect to UDG.
        fh_udg(:,1,1) = tau;                    % dfh(ne)_d(ne)
        fh_udg(:,1,5) = D_star*r.*nl(:,1);     % dfh(ne)_d(dne_dr)
        fh_udg(:,1,9) = D_star*r.*nl(:,2);     % dfh(ne)_d(dne_dz)

        fh_udg(:,2,2) = tau;                    % dfh(nn)_d(nn)
        fh_udg(:,2,6) = D_star*r.*nl(:,1);     % dfh(nn)_d(dnn_dr)
        fh_udg(:,2,10) = D_star*r.*nl(:,2);     % dfh(nn)_d(dnn_dz)

        fh_udg(:,3,3) = tau;                    % dfh(np)_d(np)
        fh_udg(:,3,7) = D_star*r.*nl(:,1);     % dfh(np)_d(dnp_dr)
        fh_udg(:,3,11) = D_star*r.*nl(:,2);     % dfh(np)_d(dnp_dz)

        fh_udg(:,4,4) = tau;
        fh_udg(:,4,8) = r.*nl(:,1);
        fh_udg(:,4,12) = r.*nl(:,2);

        % This is the jacobian of the normal flux f_hat w.r.t UHAT and does _not_ involve UDG.
        % For this problem the jacobian matrix is diagonal - each equation doesn't depend on uhat from the other equations
        fh_uh(:,1,1) = -tau;
        fh_uh(:,2,2) = -tau;
        fh_uh(:,3,3) = -tau;
        fh_uh(:,4,4) = -tau;

    case 2  % Farfield - 0 flux
        [fh,fh_udg,fh_uh] = fhat_electrondensity(nl,p,udg,uh,param,time);
        % fh = fh + ui;     We know the flux is 0 at the farfield so no need to add it in from ui

    case 3  % Outflow
        fh = zeros(ng,nch);
        fh_udg = zeros(ng,nch,nc);
        fh_uh = zeros(ng,nch,nch);

        % For now just set a dirichlet condition on all species densities. Later, add a switch for the dirichlet/neumann condition
        fh(:,1) = tau.*(0-uh(:,1));
        fh(:,2) = tau.*(0-uh(:,2));
        fh(:,3) = tau.*(0-uh(:,3));
        fh(:,4) = tau.*(0-uh(:,4));

        fh_uh(:,1,1) = -tau;
        fh_uh(:,2,2) = -tau;
        fh_uh(:,3,3) = -tau;
        fh_uh(:,4,4) = -tau;

    case 4  % Grounded surface
        fh = zeros(ng,nch);
        fh_udg = zeros(ng,nch,nc);
        fh_uh = zeros(ng,nch,nch);

        % Electrons and negative ions: no diffusive flux (absorbing boundary)

        fh(:,1) = r.*(De_star*dne_dr.*nl(:,1) + De_star*dne_dz.*nl(:,2)) + tau.*(ne-uh(:,1));
        fh(:,2) = r.*(Dn_star*dnn_dr.*nl(:,1) + Dn_star*dnn_dz.*nl(:,2)) + tau.*(nn-uh(:,2));

        fh_udg(:,1,1) = tau;                    % dfh(ne)_d(ne)
        fh_udg(:,1,5) = De_star*r.*nl(:,1);     % dfh(ne)_d(dne_dr)
        fh_udg(:,1,9) = De_star*r.*nl(:,2);     % dfh(ne)_d(dne_dz)

        fh_udg(:,2,2) = tau;                    % dfh(nn)_d(nn)
        fh_udg(:,2,6) = Dn_star*r.*nl(:,1);     % dfh(nn)_d(dnn_dr)
        fh_udg(:,2,10) = Dn_star*r.*nl(:,2);     % dfh(nn)_d(dnn_dz)

        fh_uh(:,1,1) = -tau;
        fh_uh(:,2,2) = -tau;

        % Positives: homogeneous dirichlet
        fh(:,3) = tau.*(0-uh(:,3));
        fh_uh(:,3,3) = -tau;

        % Grounded surface: homogeneous dirichlet for electrostatic
        fh(:,4) = tau.*(0-uh(:,4));
        fh_uh(:,4,4) = -tau;

    case 5  % Needle
        fh = zeros(ng,nch);
        fh_udg = zeros(ng,nch,nc);
        fh_uh = zeros(ng,nch,nch);

        normE = sqrt(Ex.^2 + Ey.^2);
        gamma = 0.1;

        % Electrons: outflow flux -- should I make this more like a neumann condition: get fhat like usual and then add the flux to that?
        % fh(:,1) = r.*(gamma.*np.*normE) + tau.*(ne-uh(:,1));
        % fh_udg(:,1,1) = tau;                    % dfh(ne)_d(ne)
        % fh_udg(:,1,3) = r.*gamma.*normE;      % The only cross-term in the boundary conditions
        % fh_uh(:,1,1) = -tau;

        % 0 flux on electrons from tip for now
        [fh_tmp,fh_udg_tmp,fh_uh_tmp] = fhat_electrondensity(nl,p,udg,uh,param,time);
        fh(:,1) = fh_tmp(:,1)+r.*(gamma.*np.*normE);
        fh_udg(:,1,:) = fh_udg_tmp(:,1,:);
        fh_uh(:,1,1) = fh_uh_tmp(:,1,1);
        fh_udg(:,1,3) = r.*gamma.*normE;      % The only cross-term in the boundary conditions
        fh_udg(:,1,8) = r.*(gamma.*np.*Ex./normE);
        fh_udg(:,1,12) = r.*(gamma.*np.*Ey./normE);

        % Negatives: homogeneous dirichlet
        fh(:,2) = tau.*(0-uh(:,2));
        fh_uh(:,2,2) = -tau;

        % Positives: no diffusive flux (absorbing boundary)
        fh(:,3) = r.*((Dp_star*dnp_dr).*nl(:,1) + (Dp_star*dnp_dz).*nl(:,2)) + tau.*(np-uh(:,3));
        fh_udg(:,3,3) = tau;                    % dfh(np)_d(np)
        fh_udg(:,3,6) = Dp_star*r.*nl(:,1);     % dfh(np)_d(dnp_dr)
        fh_udg(:,3,9) = Dp_star*r.*nl(:,2);     % dfh(np)_d(dnp_dz)
        fh_uh(:,3,3) = -tau;

        % Homogeneous dirichlet for electrostatic BY SUPERPOSITION
        fh(:,4) = tau.*(0-uh(:,4));
        fh_uh(:,4,4) = -tau;

    otherwise
        error('unknown boundary type');
end


% switch ib
%     case 1  % Dirichlet
%         fh = tau.*(ui-uh);
%         fh_udg = zeros(ng,nch,nc);
%         fh_uh = -tau;
%     case 3  % Prescribed flux - scalar
%         [fh,fh_udg,fh_uh] = fhat_electrondensity(nl,p,udg,uh,param,time);
%         fh = fh + ui;
%     case 4  % Diffusion flux only
%         r = p(:,1);
%         u = udg(:,1);
%         dne_dr = udg(:,2);
%         dne_dz = udg(:,3);
%         fh = r.*(D_star*dne_dr.*nl(:,1) + D_star*dne_dz.*nl(:,2)) + tau.*(u-uh);

%         % This is the jacobian of the normal flux f_hat w.r.t UDG. UDG is [u, uq1x, q2x], and does _not_ involve uh.
%         fh_udg = zeros(ng,nch,nc);
%         fh_udg(:,1,1) = tau;
%         fh_udg(:,1,2) = D_star*r.*nl(:,1);
%         fh_udg(:,1,3) = D_star*r.*nl(:,2);

%         % This is the jacobian of the normal flux f_hat w.r.t UHAT and does _not_ involve UDG.
%         fh_uh = zeros(ng,nch,nch);
%         fh_uh(:,1,1) = -tau; 
%     case 5 % Symmetry BC
%         r = 1.0;
%         D_star = 1.0;
%         u = udg(:,1);
%         dne_dr = udg(:,2);
%         dne_dz = udg(:,3);
%         fh = r.*(D_star*dne_dr.*nl(:,1) + D_star*dne_dz.*nl(:,2)) + tau.*(u-uh);

%         % This is the jacobian of the normal flux f_hat w.r.t UDG. UDG is [u, uq1x, q2x], and does _not_ involve uh.
%         fh_udg = zeros(ng,nch,nc);
%         fh_udg(:,1,1) = tau;
%         fh_udg(:,1,2) = D_star*r.*nl(:,1);
%         fh_udg(:,1,3) = D_star*r.*nl(:,2);

%         % This is the jacobian of the normal flux f_hat w.r.t UHAT and does _not_ involve UDG.
%         fh_uh = zeros(ng,nch,nch);
%         fh_uh(:,1,1) = -tau;
%     case 6  % Prescribed flux - function
%         [fh,fh_udg,fh_uh] = fhat_electrondensity(nl,p,udg,uh,param,time);
%         Er = udg(:,8);    % Ex is -grad(phi)
%         Ez = udg(:,12);
%         np = udg(:,2);
%         gamma = 0.001;      % Secondary emission coefficient from needle

%         normE = sqrt(Er.^2 + Ez.^2);

%         needle_flux = gamma*np*normE;
%         fh = fh + needle_flux;
%     otherwise
%         error('unknown boundary type');
% end

