function u = cghelmholtz(master,mesh,UDG,alpha)

%HDG_POSTPROCESS postprocesses the HDG solution to obtain a better solution.
%   [ustarh]=hdg_postprocess(mesh,master,uh,qh,uhath)
%
%      MASTER:       Master structure of porder
%      MESH:         Mesh structure of porder
%      MASTER1:      Master structure of porder+1
%      MESH1:        Mesh structure of porder+1
%      UH:           Approximate scalar variable
%      QH:           Approximate flux
%      USTARH:       Postprocessed scalar variable


nd    = master.nd;
npe   = master.npe;
nge   = master.nge;
ne    = size(mesh.dgnodes,3);

il = zeros(npe,npe,ne);
jl = zeros(npe,npe,ne);
for i=1:ne    
    con = mesh.t2(i,:)';    
    com = repmat(con,[1 npe]);%[con con con con con con];
    il(:,:,i) = com;
    jl(:,:,i) = com';        
end

p = mesh.p2(mesh.t2',:);
p = reshape(p,npe,ne,nd);
p = permute(p,[1 3 2]);

dshapvt   = reshape(permute(master.shapvt(:,:,2:nd+1),[1 3 2]),[nge*nd npe]);
dshapvtxi = master.shapvl(:,:,2)';
dshapvtet = master.shapvl(:,:,3)';

Ke = zeros(npe,npe,ne);
Me = zeros(npe,npe,ne);
Fe = zeros(npe,ne);
for i=1:ne
    dg  = p(:,1:nd,i);

    % compute the Jacobian matrix at Gauss points: dx/dxi
    Jg = dshapvt*dg(:,1:nd);
    Jg = reshape(Jg,[nge nd nd]);        
    [jac,~] = volgeom(Jg);
                
    dshapdx = bsxfun(@times,dshapvtxi,Jg(:,2,2)./jac)-bsxfun(@times,dshapvtet,Jg(:,1,2)./jac);
    dshapdy = bsxfun(@times,dshapvtet,Jg(:,1,1)./jac)-bsxfun(@times,dshapvtxi,Jg(:,2,1)./jac);    
    Ke(:,:,i) = dshapdx'*diag(master.gwvl.*jac)*(dshapdx)+...
                dshapdy'*diag(master.gwvl.*jac)*(dshapdy);
            
%     dshapdx  = zeros(npe,nge,nd);
%     for ii=1:nd                
%         dshapdx(:,:,ii)  = bsxfun(@times,master.shapvl(:,:,2),Xx(:,ii,1)');
%         for jj=2:nd
%             dshapdx(:,:,ii)  = dshapdx(:,:,ii) + bsxfun(@times,master.shapvl(:,:,1+jj),Xx(:,ii,jj)');
%         end            
%     end    
%     Ke(:,:,i) = dshapdx(:,:,1)*diag(master.gwvl./jac)*(dshapdx(:,:,1)');
%     for ii=2:nd
%         Ke(:,:,i) = Ke(:,:,i) + dshapdx(:,:,ii)*diag(master.gwvl./jac)*(dshapdx(:,:,ii)');
%     end    
    
    Me(:,:,i) = master.shapvl(:,:,1)*diag(master.gwvl.*jac)*master.shapvl(:,:,1)';
    ug = master.shapvl(:,:,1)'*UDG(:,1,i);
    Fe(:,i) = master.shapvl(:,:,1)*(master.gwvl.*ug.*jac);       
end

% il = zeros(npe,npe,ne);
% jl = zeros(npe,npe,ne);
% for i=1:ne    
%     con = cgelcon(:,i);    
%     com = repmat(con,[1 npe]);
%     il(:,:,i) = com;
%     jl(:,:,i) = com';        
% end

K = sparse(reshape(il,npe*npe*ne,1),reshape(jl,npe*npe*ne,1),reshape(Ke,npe*npe*ne,1));        
M = sparse(reshape(il,npe*npe*ne,1),reshape(jl,npe*npe*ne,1),reshape(Me,npe*npe*ne,1));   
F = sparse(reshape(il(:,1,:),npe*ne,1),ones(npe*ne,1),reshape(Fe,npe*ne,1)); 

u = (M + alpha*K)\F;
u = u(:);

function [jac,Xx] = volgeom(Jg)

nge = size(Jg,1); 
nd  = size(Jg,2);
switch nd
    case 1
        jac = Jg;
        Xx = -ones(nge,1);
    case 2
        jac = Jg(:,1,1).*Jg(:,2,2) - Jg(:,1,2).*Jg(:,2,1);
        Xx(:,1,1) = -Jg(:,2,2);
        Xx(:,2,1) = Jg(:,2,1);
        Xx(:,1,2) = Jg(:,1,2);
        Xx(:,2,2) = -Jg(:,1,1);
    case 3
        jac = Jg(:,1,1).*Jg(:,2,2).*Jg(:,3,3) - Jg(:,1,1).*Jg(:,3,2).*Jg(:,2,3)+ ...
              Jg(:,2,1).*Jg(:,3,2).*Jg(:,1,3) - Jg(:,2,1).*Jg(:,1,2).*Jg(:,3,3)+ ...
              Jg(:,3,1).*Jg(:,1,2).*Jg(:,2,3) - Jg(:,3,1).*Jg(:,2,2).*Jg(:,1,3);            
        Xx(:,1,1) = Jg(:,2,3).*Jg(:,3,2) - Jg(:,2,2).*Jg(:,3,3);
        Xx(:,2,1) = Jg(:,2,1).*Jg(:,3,3) - Jg(:,2,3).*Jg(:,3,1);
        Xx(:,3,1) = Jg(:,2,2).*Jg(:,3,1) - Jg(:,2,1).*Jg(:,3,2);
        Xx(:,1,2) = Jg(:,1,2).*Jg(:,3,3) - Jg(:,1,3).*Jg(:,3,2);
        Xx(:,2,2) = Jg(:,1,3).*Jg(:,3,1) - Jg(:,1,1).*Jg(:,3,3);
        Xx(:,3,2) = Jg(:,1,1).*Jg(:,3,2) - Jg(:,1,2).*Jg(:,3,1);
        Xx(:,1,3) = Jg(:,1,3).*Jg(:,2,2) - Jg(:,1,2).*Jg(:,2,3);
        Xx(:,2,3) = Jg(:,1,1).*Jg(:,2,3) - Jg(:,1,3).*Jg(:,2,1);
        Xx(:,3,3) = Jg(:,1,2).*Jg(:,2,1) - Jg(:,1,1).*Jg(:,2,2);
    otherwise
        error('Dimension is not implemented');
end




