function [UDG1,UDG2,e1,e2] = dgprojection2(master,mesh,UDG,korder,tol)
% UDG, master, mesh : porder
% UDG1 : korder

if nargin <5
    tol = 1e-6;
end

nd = mesh.nd;
dim = mesh.nd;
elemtype = mesh.elemtype;
nodetype = mesh.nodetype;
ne    = size(mesh.dgnodes,3);
nc    = size(UDG,2);
npv  = master.npv;
ngv  = master.ngv;
pgauss = 2*master.porder;

dshapvt  = reshape(permute(master.shapvt(:,:,2:dim+1),[1 3 2]),[ngv*dim npv]);

master1 = masterelement(korder,pgauss,dim,elemtype,nodetype);

npv1  = master1.npv;
UDG1 = zeros(npv1,nc,ne);
UDG2 = zeros(npv,nc,ne);
e1 = zeros(nc,ne);
e2 = zeros(nc,ne);
for i=1:ne
    dg = mesh.dgnodes(:,:,i);        
    udg = UDG(:,:,i);

    % compute the Jacobian matrix at Gauss points: dx/dxi    
    Jg = dshapvt*dg(:,1:nd);
    Jg = reshape(Jg,[ngv nd nd]);        
    jac = volgeom(Jg);        
            
    M1 = master1.shapvl(:,:,1)*diag(master.gwvl.*jac)*master1.shapvl(:,:,1)';
    C1 = master1.shapvl(:,:,1)*diag(master.gwvl.*jac)*master.shapvl(:,:,1)';
    L1 = C1*udg;
    UDG1(:,:,i) = M1\L1;        
    
    M = master.shapvl(:,:,1)*diag(master.gwvl.*jac)*master.shapvl(:,:,1)';
    C = master.shapvl(:,:,1)*diag(master.gwvl.*jac)*master1.shapvl(:,:,1)';
    L = C*UDG1(:,:,i);
    UDG2(:,:,i) = M\L;
    
    vol = sum(master.gwvl.*jac);
    edg = udg-UDG2(:,:,i);
    for j = 1:nc
        e1(j,i) = (edg(:,j)'*M*edg(:,j))/vol;        
        e2(j,i) = e1(j,i)/max(udg(:,j)'*M*udg(:,j),tol);        
    end
end

function [jac,Xx] = volgeom(Jg)

ngv = size(Jg,1); 
nd  = size(Jg,2);
switch nd
    case 1
        jac = Jg;
        Xx = -ones(ngv,1);
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




