function [fh,fh_udg,fh_uh] = fhat_cd(nl,p,udg,uh,param,time,signe)
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
nd = nq;

kappa = param{1};
c     = param{2};
tau   = param{end};

if nd==2
    r = p(:,1);
    u = udg(:,1);
    qx = udg(:,2);
    qy = udg(:,3);

    % ng x nch
    fh = r.*((c(1)*uh+kappa*qx).*nl(:,1) + (c(2)*uh+kappa*qy).*nl(:,2)) + tau.*(u-uh);

    fh_udg = zeros(ng,nch,nc);
    fh_udg(:,1,1) = tau;
    fh_udg(:,1,2) = kappa*r.*nl(:,1);
    fh_udg(:,1,3) = kappa*r.*nl(:,2);
    
    fh_uh = zeros(ng,nch,nch);
    fh_uh(:,1,1) = c(1)*r.*nl(:,1) + c(2)*r.*nl(:,2) - tau;
else
    u = udg(:,1);
    qx = udg(:,2);
    qy = udg(:,3);
    qz = udg(:,4);

    % ng x nch
    fh = kappa.*(qx.*nl(:,1)+qy.*nl(:,2)+qz.*nl(:,3)) + tau.*(u-uh);

    fh_udg = zeros(ng,nch,nc);
    fh_udg(:,1,1) = tau;
    fh_udg(:,1,2) = kappa*nl(:,1);
    fh_udg(:,1,3) = kappa*nl(:,2);    
    fh_udg(:,1,4) = kappa*nl(:,3);    
    
    fh_uh = zeros(ng,nch,nch);
    fh_uh(:,1,1) = -tau;
end


