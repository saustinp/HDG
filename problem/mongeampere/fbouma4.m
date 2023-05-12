function [fh,fh_udg,fh_uh] = fbouma4(ib,ui,nl,p,udg,uh,param,time)
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
%nq = nc-nch;

x = p(:,1);
y = p(:,2);
tau   = param{end};
%u     = udg(:,nch);
switch ib
    case 1  % Prescribed flux
        [fh,fh_udg,fh_uh] = fhatma4(nl,p,udg,uh,param,time);
        %fh = fh + ui;
    case 2 % Dirichlet
        ui = exp(0.5*(x.*x+y.*y));
        fh = tau.*(ui-uh);
        fh_udg = zeros(ng,nch,nc);
        fh_uh = -tau;     
    case 3 % Dirichlet
        ui = -sqrt(R^2 - x.^2 - y.^2);
        fh = tau.*(ui-uh);
        fh_udg = zeros(ng,nch,nc);
        fh_uh = -tau;             
    case 4 % Dirichlet
        ui = 1.0;
        fh = tau.*(ui-uh);
        fh_udg = zeros(ng,nch,nc);
        fh_uh = -tau;                     
    otherwise
        error('unknown boundary type');
end


