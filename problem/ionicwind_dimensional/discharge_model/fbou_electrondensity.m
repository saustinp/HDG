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

%%%%%%%%%       Copied from fluxgen_electrondensity.m
% Read in values from the p vector -> this normally only contains the (r,z) coordinates but we can also pass in additional fields as appended elements
r = p(:,1);
Er0 = p(:,3);       % Reading in the values of the E field under the homogeneous forcing (0 space charge) case. Keep in mind that E = -grad(phi), which is what is loaded into p(3,4)
Ez0 = p(:,4);

% Read in values from the u vector
ne = udg(:,1);
nn = udg(:,2);
np = udg(:,3);
phi = udg(:,4);
dne_dr = udg(:,5); % q is -grad(u)
dnn_dr = udg(:,6);
dnp_dr = udg(:,7);
Er_prime = udg(:,8);
dne_dz = udg(:,9);
dnn_dz = udg(:,10);
dnp_dz = udg(:,11);
Ez_prime = udg(:,12);

Er = Er_prime + Er0;
Ez = Ez_prime + Ez0;
normE = sqrt(Er.^2 + Ez.^2);

% Loading physics parameters
% The electron diffusion and mobility coefficients are functions of the reduced E field
Dn = param{3};
Dp = param{4};
E_bd = param{15};
N = param{18};
tau = param{end};

De = get_diffusion_e(normE*E_bd, N);           % Multiplying by E_bd required to convert back to dimensional units
%%%%%%%%%%%%%%%


switch ib
    case 1  % Axisymmetry boundary
        fh = zeros(ng,nch);
        fh_udg = zeros(ng,nch,nc);
        fh_uh = zeros(ng,nch,nch);

        % All equations: symmetry BC
        r = 1.0;
        D_s = 1.0;     % _s is "starred" or nondimensionalized quantity

        % The following is just copied from the fhat_electrondensity function and then the convective component taken out with all the same D_s

        % ng x nch, same size as uh
        % This is the jacobian of the normal flux f_hat w.r.t UDG. UDG is [u, uq1x, q2x], and does _not_ involve uh.
        fh(:,1) = r.*(D_s*dne_dr.*nl(:,1) + D_s*dne_dz.*nl(:,2)) + tau.*(ne-uh(:,1));
        fh(:,2) = r.*(D_s*dnn_dr.*nl(:,1) + D_s*dnn_dz.*nl(:,2)) + tau.*(nn-uh(:,2));
        fh(:,3) = r.*(D_s*dnp_dr.*nl(:,1) + D_s*dnp_dz.*nl(:,2)) + tau.*(np-uh(:,3));
        fh(:,4) = r.*(Er_prime.*nl(:,1) + Ez_prime.*nl(:,2)) + tau.*(phi-uh(:,4));

        % Jacobian of f_hat with respect to UDG.
        fh_udg(:,1,1) = tau;                    % dfh(ne)_d(ne)
        fh_udg(:,1,5) = D_s*r.*nl(:,1);     % dfh(ne)_d(dne_dr)
        fh_udg(:,1,9) = D_s*r.*nl(:,2);     % dfh(ne)_d(dne_dz)

        fh_udg(:,2,2) = tau;                    % dfh(nn)_d(nn)
        fh_udg(:,2,6) = D_s*r.*nl(:,1);     % dfh(nn)_d(dnn_dr)
        fh_udg(:,2,10) = D_s*r.*nl(:,2);     % dfh(nn)_d(dnn_dz)

        fh_udg(:,3,3) = tau;                    % dfh(np)_d(np)
        fh_udg(:,3,7) = D_s*r.*nl(:,1);     % dfh(np)_d(dnp_dr)
        fh_udg(:,3,11) = D_s*r.*nl(:,2);     % dfh(np)_d(dnp_dz)

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

        fh(:,1) = r.*(De.*dne_dr.*nl(:,1) + De.*dne_dz.*nl(:,2)) + tau.*(ne-uh(:,1));
        fh(:,2) = r.*(Dn.*dnn_dr.*nl(:,1) + Dn.*dnn_dz.*nl(:,2)) + tau.*(nn-uh(:,2));

        fh_udg(:,1,1) = tau;                    % dfh(ne)_d(ne)
        fh_udg(:,1,5) = De.*r.*nl(:,1);     % dfh(ne)_d(dne_dr)
        fh_udg(:,1,9) = De.*r.*nl(:,2);     % dfh(ne)_d(dne_dz)

        fh_udg(:,2,2) = tau;                    % dfh(nn)_d(nn)
        fh_udg(:,2,6) = Dn.*r.*nl(:,1);     % dfh(nn)_d(dnn_dr)
        fh_udg(:,2,10) = Dn.*r.*nl(:,2);     % dfh(nn)_d(dnn_dz)

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

        normE = sqrt(Er.^2 + Ez.^2);
        gamma = 0.008;

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
        fh_udg(:,1,8) = r.*(gamma.*np.*Er./normE);
        fh_udg(:,1,12) = r.*(gamma.*np.*Ez./normE);

        % Negatives: homogeneous dirichlet
        fh(:,2) = tau.*(0-uh(:,2));
        fh_uh(:,2,2) = -tau;

        % Positives: no diffusive flux (absorbing boundary)
        fh(:,3) = r.*((Dp.*dnp_dr).*nl(:,1) + (Dp.*dnp_dz).*nl(:,2)) + tau.*(np-uh(:,3));
        fh_udg(:,3,3) = tau;                    % dfh(np)_d(np)
        fh_udg(:,3,6) = Dp.*r.*nl(:,1);     % dfh(np)_d(dnp_dr)
        fh_udg(:,3,9) = Dp.*r.*nl(:,2);     % dfh(np)_d(dnp_dz)
        fh_uh(:,3,3) = -tau;

        % Homogeneous dirichlet for electrostatic BY SUPERPOSITION
        fh(:,4) = tau.*(0-uh(:,4));
        fh_uh(:,4,4) = -tau;

    otherwise
        error('unknown boundary type');
end