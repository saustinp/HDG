setapplicationpath('FM/poi');

porder = 2;
ngrid  = 33;
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
uex = exp(0.5*(x.*x+y.*y));
qex = 0*mesh.dgnodes(:,1:2,:);
vex = zeros(master.npv,mesh.nd,mesh.nd,mesh.ne);
qex(:,1,:) = x.*uex;
qex(:,2,:) = y.*uex;

app.param = [app.param 2]; 
[u,q,uhat] = hdg_poisson(master, mesh, app);

max(abs(uex(:)-u(:)))
max(abs(qex(:)-q(:)))

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
    [u,q,uhat] = hdg_poisson(master, mesh, app);
        
    eu(i,porder) = calerror(u,mesh,master,ufunc,0);    
    eq(:,i,porder) = calerror(q,mesh,master,qfunc,0);    
    
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



