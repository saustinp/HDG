function [fh,fh_udg,fh_uh] = fhat(nl,p,udg,uh,param,time,signe)
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

u     = udg(:,nch);
q     = reshape(udg(:,nch+1:nc),[ng nch nd]);

% global hx;
% global hy;

if nd ~= 1
    An = nl(:,1)*c(:,1)+nl(:,2).*c(:,2);
    
    tau = param{end}*ones(ng,1);
    %tau = 0.5*(abs(An))+kappa+1e-10;
    %tau = (abs(An))+kappa+1e-10;
%     tau = (abs(nl(:,1)) > 1.0e-2).*(kappa)*(1/hx^2) + (abs(nl(:,2)) > 1.0e-2).*(kappa)*(1/hy^2);
%     L = hx/hy
%     tau = L*(nl(:,1).*nl(:,1)) + nl(:,2).*nl(:,2);
%     tau=1000*L*ones(size(nl,1),1);
    fh = An.*uh + kappa*reshape(mapContractK(nl,q,[],2,1,3,2,1),[ng nch]) + tau.*(u-uh);

    fh_uh = An-tau;

    fh_u = tau;            
    fh_q = kappa*reshape(nl,[ng,nch,nd]);
    fh_udg = cat(3,fh_u,fh_q);
else
    An = nl(:,1)*c(:,1);
    tau = param{3};
    fh = An.*uh + kappa*reshape(mapContractK(nl,q,[],2,1,3,2,1),[ng nch]) + tau.*(u-uh);
    
    fh_uh = An-tau;
    
    fh_u = tau*ones(ng,nch,1);
    fh_q = kappa*reshape(nl,[ng,nch,nd]);
    fh_udg = cat(3,fh_u,fh_q);
end
