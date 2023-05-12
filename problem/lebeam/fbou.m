function [fh,fh_udg,fh_uh] = fbou(ib,ui,nl,p,udg,uh,param,time)
%FBOU boundary flux function

%      IC                    Boundary number
%      IB                    Boundary type
%      FUNC                  Function to compute boundary data
%      NL(N,ND)              Normal N points
%      X(N,ND)               Coordinates for N points
%      U(N,NC)               Unknown vector for N points with NC components
%      Q(N,NC,ND)            Flux vector for N points with NC components in the
%                            coordinate directions
%      P(N,1)                Pressure vector for N points with 1 component
%      M(N,NC)               Hybrid unkowns for N points with NC components
%      PARAM                 Parameter list
%      FH(N,NC):              Volume flux at N points
%      FHU(N,NC,NC):         Jacobian of the flux flux vector w.r.t. U
%      FHQ(N,NC,NC,ND):      Jacobian of the flux flux vector w.r.t. Q
%      FHP(N,NC,1):          Jacobian of the flux flux vector w.r.t. P
%      FHM(N,NC,NC):         Jacobian of the flux flux vector w.r.t. M

[ng1,nch] = size(uh);

kappa = param{1};
tau   = param{end};
tau   = tau*multiprod(ones(ng1,1),eye(nch),2,[1,2]);

uinf = repmat(ui(1,1:nch),ng1,1);

switch (ib)
    case 1 % Dirichlet bc
        fh = multiprod(tau,uinf-uh,[2,3],2);
        fh_udg = zeros(ng1,nch,3*nch);
        fh_uh = -tau;
    case 2 % homogeneous total stress bc        
        [fh,fh_udg,fh_uh] = fhatbou1(nl,p,udg,uh,param,time);    
        fh = fh + uinf;
    case 3 % homogeneous total gradient bc        
        [fh,fh_udg,fh_uh] = fhatbou2(nl,p,udg,uh,param,time);       
        fh = fh + uinf;
    case 4 % mixed Dirichlet and stress bcs
        uinf = repmat(ui(1,1),ng1,1);
        finf = repmat(ui(1,2),ng1,1);
        
        fh   = zeros(ng1,2);        
        fh_uh = zeros(ng1,nch,nch);
        fh_udg = zeros(ng1,nch,3*nch+1);
        
        fh(:,1) = tau*(uinf - uh(:,1));        
        fh_uh(:,1,1) = -tau;
        
        [gh,gh_udg,gh_uh] = fhatbou1(nl,p,udg,uh,param,time);
        fh(:,2) = gh(:,2) + finf;
        fh_uh(:,2,:)  = gh_uh(:,2,:);
        fh_udg(:,2,:) = gh_udg(:,2,:);        
    case 5 % mixed Dirichlet and stress bcs
        uinf = repmat(ui(1,1),ng1,1);
        finf = repmat(ui(1,2),ng1,1);
        
        fh   = zeros(ng1,2);        
        fh_uh = zeros(ng1,nch,nch);
        fh_udg = zeros(ng1,nch,3*nch+1);
        
        fh(:,2) = tau*(uinf - uh(:,2));        
        fh_uh(:,2,2) = -tau;
        
        [gh,gh_udg,gh_uh] = fhabout1(nl,p,udg,uh,param,time);
        fh(:,1) = gh(:,1) + finf;
        fh_uh(:,1,:)  = gh_uh(:,1,:);
        fh_udg(:,1,:) = gh_udg(:,1,:);     
    case 6 % mixed Dirichlet and gradient bcs
        uinf = repmat(ui(1,1),ng1,1);
        finf = repmat(ui(1,2),ng1,1);
        
        fh   = zeros(ng1,2);        
        fh_uh = zeros(ng1,nch,nch);
        fh_udg = zeros(ng1,nch,3*nch+1);
        
        fh(:,1) = tau*(uinf - uh(:,1));        
        fh_uh(:,1,1) = -tau;
        
        [gh,gh_udg,gh_uh] = fhatbou2(nl,p,udg,uh,param,time);
        fh(:,2) = gh(:,2) + finf;
        fh_uh(:,2,:)  = gh_uh(:,2,:);
        fh_udg(:,2,:) = gh_udg(:,2,:);        
    case 7 % mixed Dirichlet and gradient bcs
        uinf = repmat(ui(1,1),ng1,1);
        finf = repmat(ui(1,2),ng1,1);
        
        fh   = zeros(ng1,2);        
        fh_uh = zeros(ng1,nch,nch);
        fh_udg = zeros(ng1,nch,3*nch+1);
        
        fh(:,2) = tau*(uinf - uh(:,2));        
        fh_uh(:,2,2) = -tau;
        
        [gh,gh_udg,gh_uh] = fhabout2(nl,p,udg,uh,param,time);
        fh(:,1) = gh(:,1) + finf;
        fh_uh(:,1,:)  = gh_uh(:,1,:);
        fh_udg(:,1,:) = gh_udg(:,1,:);            
    otherwise
        error('unknown boundary type');
end


% total stress bc
function [fh,fh_udg,fh_uh] = fhatbou1(nl,p,udg,uh,param,time)


[ng1,nch] = size(uh);
ncs = 4;

mu  = param{1};
tau = mu;
lambda  = param{2};
cp = lambda/(mu+lambda);

fh(:,1) = mu*(2*udg(:,nch+1).*nl(:,1) + (udg(:,nch+3)+udg(:,nch+3)).*nl(:,2)) + cp*udg(:,nch+1).*nl(:,1) + tau*(udg(:,1)-uh(:,1));
fh(:,2) = mu*((udg(:,nch+2)+udg(:,nch+4)).*nl(:,1) + 2*udg(:,nch+4).*nl(:,2)) + cp*udg(:,nch+1).*nl(:,2) + tau*(udg(:,2)-uh(:,2));

if nargout > 1
    fh_u  = zeros(ng1,nch,nch);
    fh_p  = zeros(ng1,nch,1);
    fh_q  = zeros(ng1,nch,ncs);
    fh_uh = zeros(ng1,nch,nch);
    
    fh_u(:,1,1) = tau;
    fh_u(:,2,2) = tau;
    
    fh_p(:,1,1) = cp*nl(:,1);
    fh_p(:,2,1) = cp*nl(:,2);
    
    fh_q(:,1,1) = 2*mu*nl(:,1);
    fh_q(:,1,2) = mu*nl(:,2);
    fh_q(:,1,3) = mu*nl(:,2);
    fh_q(:,2,2) = mu*nl(:,1);
    fh_q(:,2,3) = mu*nl(:,1);
    fh_q(:,2,4) = 2*mu*nl(:,2);    
    
    fh_udg = cat(3,cat(3,fh_u,fh_p),fh_q);
        
    fh_uh(:,1,1) = -tau;
    fh_uh(:,2,2) = -tau;
end


% total gradient bc
function [fh,fh_udg,fh_uh] = fhatbou2(nl,p,udg,uh,param,time)

[ng1,nch] = size(uh);

kappa = param{1};

tau = kappa*multiprod(ones(ng1,1),eye(nch),2,[1,2]);

fh = kappa*(multiprod(nl(:,1),udg(:,nch+2:2*nch+1),2,2) + multiprod(nl(:,2),udg(:,2*nch+2:3*nch+1),2,2)) ...
     + [multiprod(nl(:,1),udg(:,nch+1),2,2) multiprod(nl(:,2),udg(:,nch+1),2,2)] + multiprod(tau,udg(:,1:nch)-uh,[2,3],2);
 
if nargout > 1
    fh_u = tau;
    fh_q = cat(3,multiprod(nl(:,1),kappa*eye(nch),2,[1,2]), ...
                 multiprod(nl(:,2),kappa*eye(nch),2,[1,2]));
    fh_p = [multiprod(nl(:,1),eye(1),2,[1,2]) multiprod(nl(:,2),eye(1),2,[1,2])];      
    fh_udg = cat(3,fh_u,cat(3,fh_p,fh_q));
    fh_uh = -tau;
end


