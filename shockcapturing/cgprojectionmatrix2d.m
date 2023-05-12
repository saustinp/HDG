function [M,ME] = cgprojectionmatrix2d(master, mesh, meshk)

korder = meshk.porder;

nd    = master.nd;
ne = size(mesh.t2,1);
nv = size(mesh.t2,2);
nk = size(meshk.t2,2);

p = mesh.p(mesh.t',:);
p = reshape(p,nd+1,ne,nd);
p = permute(p,[1 3 2]);

b21 = p(3,2,:)-p(1,2,:);
b22 = p(1,1,:)-p(3,1,:);
b31 = p(1,2,:)-p(2,2,:);
b32 = p(2,1,:)-p(1,1,:);
detJ = b32.*b21-b22.*b31;

il = zeros(nv,nk,ne);
jl = zeros(nv,nk,ne);
for i=1:ne    
    con = mesh.t2(i,:)';    
    com = repmat(con,[1 nk]);
    il(:,:,i) = com;
    
    con = meshk.t2(i,:)';    
    com = repmat(con,[1 nv]);
    jl(:,:,i) = com';            
end

% get shape functions and their derivatives
shapvl  = master.shapvl(:,:,1);
plocal = masternodes(korder,nd,mesh.elemtype,mesh.nodetype);
shapk = mkshape(korder,plocal,master.gpvl,mesh.elemtype);
shapvkt = shapk(:,:,1)';

ME = zeros(nv,nk,ne);
for ie = 1:ne        
    ME(:,:,ie) = shapvl*diag(master.gwvl*detJ(ie))*shapvkt;        
end

M = sparse(reshape(il,nv*nk*ne,1),reshape(jl,nv*nk*ne,1),reshape(ME,nv*nk*ne,1));        




