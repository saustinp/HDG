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

[ng,nch] = size(uh);
nc = size(udg,2);

[f,f_uhdg] = flux(p,[uh,udg(:,nch+1:end)],param,time);

tau = param{end};    
dtau_duh = 0.0;
% alpha = 100;
% s = uh.*nl(:,1) + nl(:,2);
% tau = s.*(atan(alpha*s)/pi + 0.5) - atan(alpha)/pi + 0.5;
% dtau = atan(alpha*s)/pi + (alpha*s)./(pi*(alpha^2*s.^2 + 1)) + 1/2;
% dtau_duh = dtau.*nl(:,1);

% tau = 0.5*tau;
% dtau_duh = 0.5*dtau_duh;

fh = permute(mapContractK(f,nl,2,3,1,2,[],1),[2 1]) + tau.*(udg(:,1:nch)-uh);

fh_udg = zeros(ng,nch,nc);
temp = permute(mapContractK(f_uhdg,nl,[2 4],3,1,2,[],1),[3 1 2]);
fh_udg(:,:,nch+1:end) = temp(:,:,nch+1:end);    
fh_uh = temp(:,:,1:nch);
for i=1:nch
    fh_udg(:,i,i) = tau;
    fh_uh(:,i,i) = fh_uh(:,i,i) - tau + dtau_duh.*(udg(:,i)-uh(:,i));
end    


% [ng,nc] = size(udg);
% nch = 1;
% nq = nc-nch;
% nd = nq;
% 
% kappa = param{1};
% c     = param{2};
% nu   = param{3};
% tau = param{4};
% 
% u = udg(:,1);
% qx = udg(:,2);
% qy = udg(:,3);
% 
% fh = zeros(ng,nch);
% fh_uh = zeros(ng,nch,nch);
% fh_u = zeros(ng,nch,nch);
% fh_q = zeros(ng,nch,nd);
% 
% s = qx + qy;
% alpha = 100;
% a = s.*(atan(alpha*s)/pi + 0.5) - atan(alpha)/pi + 0.5;
% dads = atan(alpha*s)/pi + (alpha*s)./(pi*(alpha^2*s.^2 + 1)) + 1/2;
% dadqx = dads;
% dadqy = dads;
% 
% if nd==2
%     fh(:,1) = (c(1)*(uh.*uh) + (nu*a+kappa).*qx).*nl(:,1) + (uh + c(2)*(uh.*uh) + (nu*a+kappa).*qy).*nl(:,2) + tau*(u-uh);
%     fh_uh(:,1,1) = (2*c(1)*uh).*nl(:,1) + (1 + 2*c(2)*uh).*nl(:,2) - tau;
% else
%     qz = udg(:,4);
%     fh(:,1) = (c(1)*(uh.*uh) + kappa*qx).*nl(:,1) + (c(2)*(uh.*uh) + kappa*qy).*nl(:,2) + (c(3)*(uh.*uh) + kappa*qz).*nl(:,3) + tau*(u-uh);
%     fh_uh(:,1,1) = (2*c(1)*uh).*nl(:,1) + (2*c(2)*uh).*nl(:,2) + (2*c(3)*uh).*nl(:,3) - tau;
%     fh_q(:,1,3) = kappa*nl(:,3);
% end
% 
% fh_u(:,1,1) = tau;
% fh_q(:,1,1) = (nu*a+kappa).*nl(:,1) + nu*dadqx.*qx.*nl(:,1) + nu*dadqx.*qy.*nl(:,2);
% fh_q(:,1,2) = (nu*a+kappa).*nl(:,2) + nu*dadqy.*qx.*nl(:,1) + nu*dadqy.*qy.*nl(:,2);
% 
% fh_udg = cat(3,fh_u,fh_q);

