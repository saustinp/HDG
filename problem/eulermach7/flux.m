function [f,f_udg] = flux(p,udg,param,time)
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
nch = 4;

gam  = param{1};
gam1 = gam - 1.0;
                                             
r    = udg(:,1);
ru   = udg(:,2);
rv   = udg(:,3);
rE   = udg(:,4);
rx   = udg(:,5);
rux  = udg(:,6);
rvx  = udg(:,7);
rEx  = udg(:,8);
ry   = udg(:,9);
ruy  = udg(:,10);
rvy  = udg(:,11);
rEy  = udg(:,12);
av = p(:,3);

r1   = 1./r;
uv   = ru.*r1;
vv   = rv.*r1;
E    = rE.*r1;
af   = 0.5*(uv.*uv+vv.*vv);
p    = gam1*(rE-r.*af);
h    = E+p.*r1;
                                        
f = zeros(ng,nch,2);
f(:,:,1) = [ru+av.*rx, ru.*uv+p+av.*rux, rv.*uv+av.*rvx,   ru.*h+av.*(rEx)];
f(:,:,2) = [rv+av.*ry, ru.*vv+av.*ruy,   rv.*vv+p+av.*rvy, rv.*h+av.*(rEy)];
% f(:,:,1) = [ru, ru.*uv+p, rv.*uv,   ru.*h];
% f(:,:,2) = [rv, ru.*vv,   rv.*vv+p, rv.*h];

if nargout>1
    f_u = zeros(ng,nch,2,nch);
    f_u(:,:,1,1) = -[zeros(ng,1), 0.5*((3-gam)*uv.*uv-gam1*vv.*vv), uv.*vv, gam*E.*uv-2*gam1*uv.*af];
    f_u(:,:,1,2) = -[-ones(ng,1), (gam-3)*uv, -vv, -gam*E+0.5*gam1*(3*uv.*uv+vv.*vv)];
    f_u(:,:,1,3) = -[zeros(ng,1), gam1*vv, -uv, gam1*uv.*vv];
    f_u(:,:,1,4) = -[zeros(ng,1), -gam1*ones(ng,1), zeros(ng,1), -gam*uv];    
    
    f_u(:,:,2,1) = -[zeros(ng,1), uv.*vv, 0.5*((3-gam)*vv.*vv-gam1*uv.*uv), gam*E.*vv-2*gam1*vv.*af];
    f_u(:,:,2,2) = -[zeros(ng,1), -vv, gam1*uv, gam1*uv.*vv];
    f_u(:,:,2,3) = -[-ones(ng,1), -uv, (gam-3)*vv,  -gam*E+0.5*gam1*(3*vv.*vv+uv.*uv) ];
    f_u(:,:,2,4) = -[zeros(ng,1), zeros(ng,1), -gam1*ones(ng,1), -gam*vv];

    %f_udg = f_u;
    
    f_q = zeros(ng,nch,2,2*nch);
    for i = 1:nch
        f_q(:,i,1,i) = av;
        f_q(:,i,2,nch+i) = av;
    end    
        
    f_udg = cat(4,f_u,f_q);
end

