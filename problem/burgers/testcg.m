porder = 4;
ngrid = 12;
hybrid = 'hdg';

% mesh and master
%mesh = mkmesh_rect(ngrid*2+1,ngrid+1,porder,0,[-1 1 0 1],1,1);
mesh = mkmesh_shockrect(porder);
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);
mesh=mkcgmesh(mesh);

i1 = find(abs(mesh.p2(:,1)-1)<1e-8);
i2 = find(abs(mesh.p2(:,1)+1)<1e-8);
i3 = find(abs(mesh.p2(:,2)-1)<1e-8);
i4 = find(abs(mesh.p2(:,2)+0)<1e-8);
mesh.ib = unique([i1;i2;i3;i4]);
mesh.in = setdiff((1:1:size(mesh.p2,1)),mesh.ib);

x1 = mesh.dgnodes(:,1,:);
x2 = mesh.dgnodes(:,2,:);
s = (2*pi*pi)*sin(pi*x1).*sin(pi*x2);
u = cgpoisson(mesh,master,s,[1 0]);

[ne, npe] = size(mesh.t2);
p = u(mesh.t2',:);
p = reshape(p,npe,1,ne);

x1 = mesh.dgnodes(:,1,:);
x2 = mesh.dgnodes(:,2,:);
ue = sin(pi*x1).*sin(pi*x2);

figure(2); clf;scaplot(mesh,p,[],2,1); axis equal; axis tight;



% mesh.ib = [];
% mesh.in = 1:size(mesh.p2,1);
% u = cgpoisson(mesh,master,s,[0.0001 1.0]);
% [ne, npe] = size(mesh.t2);
% u = u(mesh.t2');
% u = reshape(u,npe,1,ne);
% figure(2); clf;scaplot(mesh,u,[0 20],2,0); axis equal; axis tight;


