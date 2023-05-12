function [f,f_udg] = fluxma4(p,udg,param,time)
%FLUX Volume flux function
%   [f,fu,fq] = flux(p,u,q,param)
%
%      P(N,ND)              Coordinates for N points
%      U(N,NC)              Unknown vector for N points with NC components
%      Q(N,NC,ND)           Flux vector for N points with NC components in the
%                           coordinate directions
%      PARAM                Parameter list
%      F(N,NC,ND):          Volume flux at N points
%      FU(N,NC,ND,NC):      Jacobian of the flux flux vector w.r.t. U
%      FQ(N,NC,ND,NC,ND):   Jacobian of the flux flux vector w.r.t. Q

[ng,nc] = size(udg);
nch = 1;
nq = nc-nch;
nd = nq;

v11 = p(:,3);
v12 = p(:,4);
v21 = p(:,5);
v22 = p(:,6);

if nd==2
    %u = udg(:,1);
    qx = udg(:,2);
    qy = udg(:,3);

    f = zeros(ng,nch,nd);
    f(:,1,1) = v22.*qx - v21.*qy;
    f(:,1,2) = v11.*qy - v12.*qx;

    f_udg = zeros(ng,nch,nd,nc);
    f_udg(:,1,1,2) =  v22;
    f_udg(:,1,1,3) = -v21;
    f_udg(:,1,2,2) = -v12;
    f_udg(:,1,2,3) =  v11;
else
  error("nd=3 is not supported");
end

