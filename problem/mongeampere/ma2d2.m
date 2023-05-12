setapplicationpath('FM/poi');

porder = 2;
ngrid  = 17;
elemtype = 0;
nodetype = 1;
hybrid = 'hdg';

R = sqrt(2)+0.01;
app.tau = 1;
app.bcm = [1;1;1;1]*3;
app.param = R;
app.fbou = 'ubouma2';
app.source = 'sourcema2';

mesh   = mkmesh_square(ngrid,ngrid,porder,0,1,1,elemtype,nodetype);
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);
 
x = (mesh.dgnodes(:,1,:));
y = (mesh.dgnodes(:,2,:));

% f= exp(x^2 + y^2)*(x^2 + y^2 + 1)
uex = -sqrt(R^2 - x.^2 - y.^2);
qex = 0*mesh.dgnodes(:,1:2,:);
vex = zeros(master.npv,mesh.nd,mesh.nd,mesh.ne);
qex(:,1,:) = -x./uex;
qex(:,2,:) = -y./uex;
vex(:,1,1,:) = (R^2 - y.^2)./(R^2 - x.^2 - y.^2).^(3/2);
vex(:,2,2,:) = (R^2 - x.^2)./(R^2 - x.^2 - y.^2).^(3/2);
vex(:,1,2,:) = (x.*y)./(R.^2 - x.^2 - y.^2).^(3/2);
vex(:,2,1,:) = (x.*y)./(R.^2 - x.^2 - y.^2).^(3/2);

ufunc = @(x,t) -sqrt(R^2 - x(:,1).^2 - x(:,2).^2);
qfunc = @(x,t) x./sqrt(R^2 - x(:,1).^2 - x(:,2).^2);
vfunc = @(x,t) [(R^2 - x(:,2).^2) x(:,1).*x(:,2) x(:,1).*x(:,2) (R^2 - x(:,1).^2)]./(R^2 - x(:,1).^2 - x(:,2).^2).^(3/2);

u = (x.^2+y.^2)/2;
uhat = inituhat(master,mesh.elcon,u,1);
uhat = uhat(:);
[u,q,uhat,v] = hdg_ma3(master, mesh, app, u, uhat);

max(abs(uex(:)-u(:)))
max(abs(qex(:)-q(:)))
max(abs(vex(:)-v(:)))

% figure(1); clf; scaplot(mesh,uex(:,1,:),[],0,1); axis equal; axis tight; colormap jet;
% figure(2); clf; scaplot(mesh,u(:,1,:),[],0,1); axis equal; axis tight; colormap jet;

iter = zeros(5,3);
eu = zeros(5,3);
eustar = zeros(5,3);
eq = zeros(2,5,3);
ev = zeros(4,5,3);
nn = [4 8 16 32 64 128]+1;
for porder = 1:3
  for i = 1:6
    ngrid = nn(i);
    mesh   = mkmesh_square(ngrid,ngrid,porder,0,1,1,elemtype,nodetype);
    master = mkmaster(mesh,2*porder);
    [master,mesh] = preprocess(master,mesh,hybrid);
    x = (mesh.dgnodes(:,1,:));
    y = (mesh.dgnodes(:,2,:));
    u = (x.^2+y.^2)/2;
    uhat = inituhat(master,mesh.elcon,u,1);
    uhat = uhat(:);
    [u,q,uhat,v,iter(i,porder)] = hdg_ma2(master, mesh, app, u, uhat);    
    
    eu(i,porder) = calerror(u,mesh,master,ufunc,0);    
    eq(:,i,porder) = calerror(q,mesh,master,qfunc,0);    
    v1 = 0*q;
    v1(:,1,:) = v(:,1,1,:);
    v1(:,2,:) = v(:,2,1,:);
    v1(:,3,:) = v(:,1,2,:);
    v1(:,4,:) = v(:,2,2,:);
    ev(:,i,porder) = calerror(v1,mesh,master,vfunc,0);
    
    mesh1   = mkmesh_square(ngrid,ngrid,porder+1,0,1,1,elemtype,nodetype);
    master1 = mkmaster(mesh1,2*(porder+1));
    [master1,mesh1] = preprocess(master1,mesh1,hybrid);
    UDG = 0*q;
    UDG(:,1,:) = u(:,1,:);
    UDG(:,2,:) = -q(:,1,:);
    UDG(:,3,:) = -q(:,2,:);
    ustar = postprocessnd(master,mesh,master1,mesh1,UDG);
    eustar(i,porder) = calerror(ustar,mesh1,master1,ufunc,0);    
    
    [i eu(i,porder) eustar(i,porder)]
  end
end

fu=eu; cu=log(fu(1:end-1,:)./fu(2:end,:))/log(2);
fu=eustar; custar=log(fu(1:end-1,:)./fu(2:end,:))/log(2);

fq=sqrt(squeeze(eq(1,:,:).^2+eq(2,:,:).^2));
cq=log(fq(1:end-1,:)./fq(2:end,:))/log(2);

fv=sqrt(squeeze(ev(1,:,:).^2+ev(2,:,:).^2+ev(3,:,:).^2+ev(4,:,:).^2));
cv=log(fv(1:end-1,:)./fv(2:end,:))/log(2);

return;

