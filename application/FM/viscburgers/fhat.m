function [fh,fh_udg,fh_uh] = fhat(nl,p,udg,uh,param,time)
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

u = udg(:,1);
qx = udg(:,2);
qy = udg(:,3);

fh = zeros(ng,nch);
fh_uh = zeros(ng,nch,nch);
fh_u = zeros(ng,nch,nch);
fh_q = zeros(ng,nch,nd);

if nd==2
    fh(:,1) = (c(1)*(uh.*uh) + kappa*qx).*nl(:,1) + (c(2)*(uh.*uh) + kappa*qy).*nl(:,2) + tau*(u-uh);
    fh_uh(:,1,1) = (2*c(1)*uh).*nl(:,1) + (2*c(2)*uh).*nl(:,2) - tau;
else
    qz = udg(:,4);
    fh(:,1) = (c(1)*(uh.*uh) + kappa*qx).*nl(:,1) + (c(2)*(uh.*uh) + kappa*qy).*nl(:,2) + (c(3)*(uh.*uh) + kappa*qz).*nl(:,3) + tau*(u-uh);
    fh_uh(:,1,1) = (2*c(1)*uh).*nl(:,1) + (2*c(2)*uh).*nl(:,2) + (2*c(3)*uh).*nl(:,3) - tau;
    fh_q(:,1,3) = kappa*nl(:,3);
end

fh_u(:,1,1) = tau;
fh_q(:,1,1) = kappa*nl(:,1);
fh_q(:,1,2) = kappa*nl(:,2);

fh_udg = cat(3,fh_u,fh_q);

