setapplicationpath('FM/poi');

porder = 3;
ngrid  = 51;
elemtype = 0;
nodetype = 1;
hybrid = 'hdg';

R = sqrt(2)+0.01/2;
app.tau = 1;
app.bcm = [1;1;1;1]*4;
app.param = R;
app.fbou = 'ubouma3';
app.source = 'sourcema3';

mesh   = mkmesh_rect2(ngrid,ngrid,porder,0,[-1 1 -1 1],elemtype,nodetype);
%mesh = mkmesh_rect(m,n,porder,parity,xrect,elemtype,nodetype)
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);
 
x = (mesh.dgnodes(:,1,:));
y = (mesh.dgnodes(:,2,:));
u = (x.^2+y.^2)/2 + 0.3;
uhat = inituhat(master,mesh.elcon,u,1);
uhat = uhat(:);
[u,q,uhat,v] = hdg_ma2(master, mesh, app, u, uhat);
figure(1); clf; scaplot(mesh,u(:,1,:),[],0,1); axis equal; axis tight; colormap jet;

return;

