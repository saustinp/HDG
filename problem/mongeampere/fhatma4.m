function [fh,fh_udg,fh_uh] = fhatma4(nl,p,udg,uh,param,time,signe)
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

v11 = p(:,3);
v12 = p(:,4);
v21 = p(:,5);
v22 = p(:,6);

tau   = param{end};
if nd==2
    u = udg(:,1);
    qx = udg(:,2);
    qy = udg(:,3);

    % ng x nch
    fx = v22.*qx - v21.*qy;
    fy = v11.*qy - v12.*qx;
    
    fh = (fx.*nl(:,1) + fy.*nl(:,2)) + tau.*(u-uh);

    fh_udg = zeros(ng,nch,nc);
    fh_udg(:,1,1) = tau;
    fh_udg(:,1,2) = v22.*nl(:,1) - v12.*nl(:,2);
    fh_udg(:,1,3) = v11.*nl(:,2) - v21.*nl(:,1);
else
    error("nd=3 is not supported");
end

fh_uh = zeros(ng,nch,nch);
fh_uh(:,1,1) = -tau;

