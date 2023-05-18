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
nch = 1;
nq = nc-nch;

% Physics parameters
% r0 = param{1};
% z0 = param{2};
% s0 = param{3};
% Nmax = param{4};
% e = param{5};
% epsilon0 = param{6};
% Ua = param{7};
% gamma = param{8};
% E_bd = param{9};
% r_tip = param{10};
% n_ref = param{11};
% N = param{12};
% mue_ref = param{13};
D_star = param{14};
tau = param{end};
u = udg(:,nch);     % Why is u defined multiple times?

switch ib
    case 1  % Dirichlet
        fh = tau.*(ui-uh);
        fh_udg = zeros(ng,nch,nc);
        fh_uh = -tau;
    case 2  % "Extrapolate  m = u 
        fh = tau.*(u-uh);
        fh_u = tau*ones(ng,nch,nch);
        fh_q = zeros(ng,nch,nq);
        fh_udg = cat(3,fh_u,fh_q);
        fh_uh = -tau*ones(ng,nch,nch);
    case 3  % Prescribed flux
        [fh,fh_udg,fh_uh] = fhat_electrondensity(nl,p,udg,uh,param,time);
        fh = fh + ui;
    case 4  % Diffusion flux only
        r = p(:,1);
        u = udg(:,1);
        dne_dr = udg(:,2);
        dne_dz = udg(:,3);
        fh = r.*(D_star*dne_dr.*nl(:,1) + D_star*dne_dz.*nl(:,2)) + tau.*(u-uh);

        fh_udg = zeros(ng,nch,nc);
        fh_udg(:,1,1) = tau;
        fh_udg(:,1,2) = D_star*r.*nl(:,1);
        fh_udg(:,1,3) = D_star*r.*nl(:,2);

        fh_uh = zeros(ng,nch,nch);
        fh_uh(:,1,1) = -tau; 
    case 5 % Symmetry boundary condition
        r = 1.0;
        D_star = 1.0;
        u = udg(:,1);
        dne_dr = udg(:,2);
        dne_dz = udg(:,3);
        fh = r.*(D_star*dne_dr.*nl(:,1) + D_star*dne_dz.*nl(:,2)) + tau.*(u-uh);
        fh_udg = zeros(ng,nch,nc);
        fh_udg(:,1,1) = tau;
        fh_udg(:,1,2) = D_star*r.*nl(:,1);
        fh_udg(:,1,3) = D_star*r.*nl(:,2);
        fh_uh = zeros(ng,nch,nch);
        fh_uh(:,1,1) = -tau;
    % case 6  % Prescribed flux according to function or other properties
    %     [fh,fh_udg,fh_uh] = fhat_electrondensity(nl,p,udg,uh,param,time);
    %     r = p(:,1);
    %     Ex = p(:,3);    % Ex is -grad(phi)
    %     Ey = p(:,4);
    %     np = udg(:,??);     % Need to figure out where in the vector the np goes
    %     gamma = 0.001;      % Secondary emission coefficient from needle

    %     normE = sqrt(Ex.^2 + Ey.^2);

    %     needle_flux = gamma*np*normE;

    %     fh = fh + needle_flux;
    
    % case 6 % Dirichlet - custon function
    %     x = p(:,1);
    %     y = p(:,2);
    %     ui = sin(0.5*pi*x).*sin(0.5*pi*y);       
    %     fh = tau.*(ui-uh);
    %     fh_udg = zeros(ng,nch,nc);
    %     fh_uh = -tau;     
    otherwise
        error('unknown boundary type');
end


