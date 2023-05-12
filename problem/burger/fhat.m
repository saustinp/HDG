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

[fh,fh_udg,fh_uh] = fhat0(nl,p,udg,uh,param,time);

end


function [fh,fh_udg,fh_uh] = fhat0(nl,p,udg,uh,param,time)

[ng,nch] = size(uh);
nc = size(udg,2);

[f,f_uhdg] = flux(p,[uh,udg(:,nch+1:end)],param,time);

tau = 1;
dtau_du = 0;
dtau_duh = 0;
 
fh = permute(mapContractK(f,nl,2,3,1,2,[],1),[2 1]) + tau.*(udg(:,1:nch)-uh);

fh_udg = zeros(ng,nch,nc);
temp = permute(mapContractK(f_uhdg,nl,[2 4],3,1,2,[],1),[3 1 2]);
fh_udg(:,:,nch+1:end) = temp(:,:,nch+1:end);    
fh_uh = temp(:,:,1:nch);
for i=1:nch
    fh_udg(:,i,i) = fh_udg(:,i,i) + tau + dtau_du.*(udg(:,i)-uh(:,i));
    fh_uh(:,i,i) = fh_uh(:,i,i) - tau + dtau_duh.*(udg(:,i)-uh(:,i));
end    

end

function [fh,fh_udg,fh_uh] = fhat1(nl,p,udg,uh,param,time)

[ng,nch] = size(uh);
nc = size(udg,2);

[f,f_uhdg] = flux(p,[uh,udg(:,nch+1:end)],param,time);

%s = uh(:,1).*nl(:,1) + nl(:,2);
alpha = 1000;
s = uh.*nl(:,1) + nl(:,2);
tau = s.*(atan(alpha*s)/pi + 0.5) - atan(alpha)/pi + 0.5;
dtau = atan(alpha*s)/pi + (alpha*s)./(pi*(alpha^2*s.^2 + 1)) + 1/2;
dtau_duh = dtau.*nl(:,1);
dtau_du = 0;

% tau = 0.5*tanh(alpha*s)+0.5;
% dtau = -(alpha*(tanh(alpha*s).^2 - 1))/2;
% dtau_du = 0;
% dtau_duh = dtau.*nl(:,1);
% 
fh = permute(mapContractK(f,nl,2,3,1,2,[],1),[2 1]) + tau.*(udg(:,1:nch)-uh);

fh_udg = zeros(ng,nch,nc);
temp = permute(mapContractK(f_uhdg,nl,[2 4],3,1,2,[],1),[3 1 2]);
fh_udg(:,:,nch+1:end) = temp(:,:,nch+1:end);    
fh_uh = temp(:,:,1:nch);
for i=1:nch
    fh_udg(:,i,i) = fh_udg(:,i,i) + tau + dtau_du.*(udg(:,i)-uh(:,i));
    fh_uh(:,i,i) = fh_uh(:,i,i) - tau + dtau_duh.*(udg(:,i)-uh(:,i));
end    

end

function [fh,fh_udg,fh_uh] = fhat2(nl,p,udg,uh,param,time)

[ng,nch] = size(uh);
nc = size(udg,2);

u = udg(:,1:nch);

[f1,f1_uhdg] = flux(p,[u,udg(:,nch+1:end)],param,time);
[f2,f2_uhdg] = flux(p,[2*uh-u,udg(:,nch+1:end)],param,time);
f = 0.5*(f1+f2);

s = uh.*nl(:,1) + nl(:,2);
tau = 2*abs(s);
dtau_duh = 2*sign(s)*nl(:,1);
epslm = 0.005;
if epslm>0
    rlam = (tau.*tau./(epslm)+epslm);
    ic = tau < epslm;
    tau = ic.*rlam + (1-ic).*tau;
    dtau_duh = ic.*(2.*tau.nl(:,1)/(epslm)) + (1-ic).*2*sign(s)*nl(:,1);
end
dtau_du = 0;
%dtau_duh = 0*nl(:,1);

% idx = s<0;
% dtau_duh(idx) = -nl(idx,1);

% alpha = 1000;
% s = uh.*nl(:,1) + nl(:,2);
% tau = s.*(atan(alpha*s)/pi + 0.5) - atan(alpha)/pi + 0.5;
% dtau = atan(alpha*s)/pi + (alpha*s)./(pi*(alpha^2*s.^2 + 1)) + 1/2;
% dtau_duh = dtau.*nl(:,1);
% dtau_du = 0;

% tau = 0.5*tanh(alpha*s)+0.5;
% dtau = -(alpha*(tanh(alpha*s).^2 - 1))/2;
% dtau_du = 0;
% dtau_duh = dtau.*nl(:,1);

% tau = 1;
% dtau_duh = 0;
% dtau_du = 0;
% 
% s = 0.5.*nl(:,1) + 1.0.*nl(:,2);
% alpha = 1e4;
% tau = s.*(atan(alpha*s)/pi + 0.5) - atan(alpha)/pi + 0.5;
% tau = max(1e-4,s);
% % tau = abs(s)+1e-4;

temp1 = permute(mapContractK(f1_uhdg,nl,[2 4],3,1,2,[],1),[3 1 2]);
temp2 = permute(mapContractK(f2_uhdg,nl,[2 4],3,1,2,[],1),[3 1 2]);

% t = 0.5*(temp1(:,:,nch+1:end) + temp2(:,:,nch+1:end));
% max(abs(t(:)))

fh = permute(mapContractK(f,nl,2,3,1,2,[],1),[2 1]) + tau.*(udg(:,1:nch)-uh);
fh_udg = zeros(ng,nch,nc);
fh_udg(:,:,nch+1:end) = 0.5*(temp1(:,:,nch+1:end) + temp2(:,:,nch+1:end));    
fh_uh = temp2(:,:,1:nch);
fh_udg(:,:,1:nch) = 0.5*(temp1(:,:,1:nch)-temp2(:,:,1:nch));
for i=1:nch
    fh_udg(:,i,i) = fh_udg(:,i,i) + tau + dtau_du.*(udg(:,i)-uh(:,i));
    fh_uh(:,i,i) = fh_uh(:,i,i) - tau + dtau_duh.*(udg(:,i)-uh(:,i));
end    

end


% [ng,nch] = size(uh);
% nc = size(udg,2);
% 
% [f,f_udg] = flux(p,udg(:,1:end),param,time);
% %[f,f_udg] = flux(p,[uh,udg(:,nch+1:end)],param,time);
% 
% %tau = param{end};    
% alpha = 100;
% s = uh.*nl(:,1) + nl(:,2);
% tau = s.*(atan(alpha*s)/pi + 0.5) - atan(alpha)/pi + 0.5;
% dtau = atan(alpha*s)/pi + (alpha*s)./(pi*(alpha^2*s.^2 + 1)) + 1/2;
% dtau_duh = dtau.*nl(:,1);
% 
% % tau = 0.5*tau;
% % dtau_duh = 0.5*dtau_duh;
% 
% % tau = ones(ng,1);
% % dtau_duh = zeros(ng,1);
% % idx = find(s>0);
% % tau(idx) = 0;
% 
% fh = permute(mapContractK(f,nl,2,3,1,2,[],1),[2 1]) + tau.*(udg(:,1:nch)-uh);
% 
% fh_udg = zeros(ng,nch,nc);
% temp = permute(mapContractK(f_udg,nl,[2 4],3,1,2,[],1),[3 1 2]);
% fh_udg(:,:,1:end) = temp(:,:,1:end);    
% fh_uh = 0*temp(:,:,1:nch);
% for i=1:nch
%     fh_udg(:,i,i) = fh_udg(:,i,i) + tau;
%     fh_uh(:,i,i) = fh_uh(:,i,i) - tau + dtau_duh.*(udg(:,i)-uh(:,i));
% end    
% 
% 
% [f,f_uhdg] = flux(p,[uh,udg(:,nch+1:end)],param,time);
% 
% %tau = param{end};    
% alpha = 1000;
% s = uh.*nl(:,1) + nl(:,2);
% tau = s.*(atan(alpha*s)/pi + 0.5) - atan(alpha)/pi + 0.5;
% dtau = atan(alpha*s)/pi + (alpha*s)./(pi*(alpha^2*s.^2 + 1)) + 1/2;
% dtau_duh = dtau.*nl(:,1);
% 
% % tau = 0.5*tau;
% % dtau_duh = 0.5*dtau_duh;
% 
% % tau = ones(ng,1);
% % dtau_duh = zeros(ng,1);
% % idx = find(s>0);
% % tau(idx) = 0;
% 
% fh = permute(mapContractK(f,nl,2,3,1,2,[],1),[2 1]) + tau.*(udg(:,1:nch)-uh);
% 
% fh_udg = zeros(ng,nch,nc);
% temp = permute(mapContractK(f_uhdg,nl,[2 4],3,1,2,[],1),[3 1 2]);
% fh_udg(:,:,nch+1:end) = temp(:,:,nch+1:end);    
% fh_uh = temp(:,:,1:nch);
% for i=1:nch
%     fh_udg(:,i,i) = tau;
%     fh_uh(:,i,i) = fh_uh(:,i,i) - tau + dtau_duh.*(udg(:,i)-uh(:,i));
% end    
% 
