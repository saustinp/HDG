setapplicationpath('FM/poi');

porder = 1;
ngrid  = 5;
elemtype = 1;
nodetype = 1;
hybrid = 'hdg';

app.tau = 1;
app.bcm = [1;1;1;1]*0;
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
[u,q,uhat,v,iter] = hdg_ma(master, mesh, app, u, uhat);

max(abs(uex(:)-u(:)))
max(abs(qex(:)-q(:)))

app.param = [app.param 3]; 
vdg = 0*q;
vdg(:,1,:) = v(:,1,1,:);
vdg(:,2,:) = v(:,2,2,:);
[u,q,uhat] = hdg_poisson(master, mesh, app, vdg);

max(abs(uex(:)-u(:)))
max(abs(qex(:)-q(:)))
max(abs(vex(:)-v(:)))

for i = 1:2
  for j=1:2
    tm = vex(:,i,j,:)-v(:,i,j,:);
    max(abs(tm(:)))
  end
end

figure(1); clf; scaplot(mesh,uex(:,1,:),[],0,1); axis equal; axis tight; colormap jet;
figure(2); clf; scaplot(mesh,u(:,1,:),[],0,1); axis equal; axis tight; colormap jet;

ufunc = @(x,t) exp(0.5*(x(:,1).*x(:,1)+x(:,2).*x(:,2)));
eu = calerror(u,mesh,master,ufunc,0);

qfunc = @(x,t) x.*(exp(0.5*(x(:,1).*x(:,1)+x(:,2).*x(:,2))));
eq = calerror(q,mesh,master,qfunc,0);

vfunc = @(x,t) [1+x(:,1).*x(:,1) x(:,1).*x(:,2) x(:,1).*x(:,2) 1+x(:,2).*x(:,2)].*(exp(0.5*(x(:,1).*x(:,1)+x(:,2).*x(:,2))));
v1 = 0*q;
v1(:,1,:) = v(:,1,1,:);
v1(:,2,:) = v(:,2,1,:);
v1(:,3,:) = v(:,1,2,:);
v1(:,4,:) = v(:,2,2,:);
ev = calerror(v1,mesh,master,vfunc,0);


mesh1   = mkmesh_square(ngrid,ngrid,porder+1,0,1,1,elemtype,nodetype);
master1 = mkmaster(mesh1,2*(porder+1));
[master1,mesh1] = preprocess(master1,mesh1,hybrid);

UDG = 0*q;
UDG(:,1,:) = u(:,1,:);
UDG(:,2,:) = -q(:,1,:);
UDG(:,3,:) = -q(:,2,:);
ustar = postprocessnd(master,mesh,master1,mesh1,UDG);

UDG(:,1,:) = q(:,1,:);
UDG(:,2,:) = -v(:,1,1,:);
UDG(:,3,:) = -v(:,2,1,:);
qstar = mesh1.dgnodes;
qstar(:,1,:) = postprocessnd(master,mesh,master1,mesh1,UDG);

UDG(:,1,:) = q(:,2,:);
UDG(:,2,:) = -v(:,1,2,:);
UDG(:,3,:) = -v(:,2,2,:);
qstar(:,2,:) = postprocessnd(master,mesh,master1,mesh1,UDG);

x = (mesh1.dgnodes(:,1,:));
y = (mesh1.dgnodes(:,2,:));
ue = exp(0.5*(x.*x+y.*y));
qe = mesh1.dgnodes;
qe(:,1,:) = x.*ue;
qe(:,2,:) = y.*ue;

max(abs(ue(:)-ustar(:)))
max(abs(qe(:)-qstar(:)))

iter = zeros(5,3);
eu = zeros(5,3);
eustar = zeros(5,3);
eq = zeros(2,5,3);
ev = zeros(4,5,3);
nn = [4 8 16 32 64]+1;
for porder = 1:3
  for i = 1:5
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


% 
% fq=sqrt(squeeze(err(2,:,:).^2+err(3,:,:).^2));
% fx=errstar;
% 
% cu=log(fu(1:end-1,:)./fu(2:end,:))/log(2);
% cq=log(fq(1:end-1,:)./fq(2:end,:))/log(2);
% cx=log(fx(1:end-1,:)./fx(2:end,:))/log(2);


return;

