porder = 4;
ngrid  = 21;
elemtype = 1;
nodetype = 1;
hybrid = 'hdg';

mesh   = mkmesh_square(ngrid,ngrid,porder,0,1,1,elemtype,nodetype);
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);
mesh=mkcgmesh(mesh); 
i1 = find(abs(mesh.p2(:,1)-1)<1e-8);
i2 = find(abs(mesh.p2(:,1)+0)<1e-8);
i3 = find(abs(mesh.p2(:,2)-1)<1e-8);
i4 = find(abs(mesh.p2(:,2)+0)<1e-8);
mesh.ib = unique([i1;i2;i3;i4]);
mesh.in = setdiff((1:1:size(mesh.p2,1)),mesh.ib);

x1 = mesh.p2(:,1);
x2 = mesh.p2(:,2);
s = (2*pi*pi)*sin(pi*x1).*sin(pi*x2);
u = cgpoisson(mesh,master,s,[1 0]);

[ne, npe] = size(mesh.t2);
p = u(mesh.t2',:);
p = reshape(p,npe,1,ne);

x1 = mesh.dgnodes(:,1,:);
x2 = mesh.dgnodes(:,2,:);
ue = sin(pi*x1).*sin(pi*x2);

figure(2); clf;scaplot(mesh,p-ue,[],2,0); axis equal; axis tight;


korder = 3;
meshk  = mkmesh_square(ngrid,ngrid,korder,0,1,1,elemtype,nodetype);
meshk = mkcgmesh(meshk); 
ent2elem = mkent2elem(meshk.t2');
nent = max(meshk.t2(:));
nelem = zeros(nent,1);
for i = 1:nent % for CG node i
    nelem(i) = ent2elem(i+1)-ent2elem(i); % number of elements connected to CG node i    
end

M1 = dgprojectionmatrix(meshk,korder);
M2 = dgprojectionmatrix(mesh, korder);
u1 = zeros(size(M1,1),1,meshk.ne);
ucg = zeros(nent,1);
for i = 1:mesh.ne
    u1(:,:,i) = M1(:,:,i)\M2(:,:,i)*p(:,1,i);
    ucg(meshk.t2(i,:),:) = ucg(meshk.t2(i,:),:) + u1(:,:,i); 
end
ucg = ucg./nelem;
u2 = ucg(meshk.t2');

figure(2); clf;scaplot(meshk,u1,[],2,0); axis equal; axis tight;
figure(3); clf;scaplot(meshk,u2-squeeze(u1),[],2,0); axis equal; axis tight;

% meshk   = mkmesh_square(ngrid,ngrid,korder,0,1,1,elemtype,nodetype);
% meshk =mkcgmesh(meshk); 
% [M,ME] = cgprojectionmatrix2d(master, mesh, meshk);

% dim = 2;
% porder2=5; 
% 
% plocal = masternodes(porder,dim,elemtype,nodetype);
% A = VandermondeMatrix(porder, plocal, elemtype);
% %A = inv(A);
% 
% plocal2 = masternodes(porder2,dim,elemtype,nodetype);
% A2 = VandermondeMatrix(porder2, plocal2, elemtype);
% 
% mesh2   = mkmesh_square(ngrid,ngrid,porder2,0,1,1,elemtype,nodetype);
% master2 = mkmaster(mesh2,2*porder2);
% [~,mesh2] = preprocess(master2,mesh2,hybrid);
% 
% [~,~,rowent2elem,~,cgent2dgent] = mkcgent2dgent(mesh2.dgnodes(:,1:2,:),1e-8);
% 
% [u2, udg2] = cgprojection(u, A, A2, mesh.t2', cgent2dgent, rowent2elem); 
% 
% figure(3); clf;scaplot(mesh2,udg2,[],2,1); axis equal; axis tight;


