setapplicationpath('FM/poi');

porder = 5;
ngrid  = 5;
elemtype = 0;
nodetype = 1;
hybrid = 'hdg';

app.bcm = [1;1;1;1]*2;
app.param = [kappa,tau];
app.fbou = 'ubouma';
app.source = 'sourcema';

mesh   = mkmesh_square(ngrid,ngrid,porder,0,1,1,elemtype,nodetype);
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);
 
x = (mesh.dgnodes(:,1,:));
y = (mesh.dgnodes(:,2,:));

% f= exp(x^2 + y^2)*(x^2 + y^2 + 1)
uex = exp(0.5*(x.*x+y.*y));
qex = 0*mesh.dgnodes(:,1:2,:);
vex = zeros(master.npv,mesh.nd,mesh.nd,mesh.ne);
qex(:,1,:) = x.*uex;
qex(:,2,:) = y.*uex;
vex(:,1,1,:) = uex + x.^2.*uex;
vex(:,2,2,:) = uex + y.^2.*uex;
vex(:,1,2,:) = x.*y.*uex;
vex(:,2,1,:) = x.*y.*uex;

u = (x.^2+y.^2)/2;
uhat = inituhat(master,mesh.elcon,u,1);
uhat = uhat(:);
[u,q,uhat,v] = hdg_ma(master, mesh, app, u, uhat);

max(abs(uex(:)-u(:)))
max(abs(qex(:)-q(:)))
max(abs(vex(:)-v(:)))

figure(1); clf; scaplot(mesh,uex(:,1,:),[],0,1); axis equal; axis tight; colormap jet;
figure(2); clf; scaplot(mesh,u(:,1,:),[],0,1); axis equal; axis tight; colormap jet;

return;

