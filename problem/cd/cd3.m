setapplicationpath('FM/condiff');

porder = 5;
ngrid  = 17;
hybrid = 'hdg';

kappa = 1e-3;
c = [1,0]; 
tau = 2;

app.adjoint = 0;
app.denseblock = 0;
app.hybrid = 'hdg';
app.localsolve=1;
app.arg = {kappa,c,tau};
app.bcm = [1;4;1;3];
app.bcs = [1;0;0;0]; 

% for adjoint problem
app.bcd = [4;1;1;1];
app.bcv = [1;0;0;0]; 

app.tdep = false;
app.wave = false;
app.alag = false;
app.flg_q = 1;
app.flg_p = 0;
app.flg_g = 0;

app.fc_q = 1;
app.fc_u = 0;
app.fc_p = 0;

app.np = 2;
app.nd = 2;
app.nch  = 1;                       % Number of componets of UH
app.nc   = app.nch*(app.nd+1);    % Number of componeents of UDG
app.ncu = 1;

app.time = [];
app.dtfc = [];
app.alpha = [];

%mesh   = mkmesh_square(ngrid,ngrid,porder);
%mesh = mkmesh_rect(ngrid,ngrid,porder,0,[-1.5 1.5 0 1],0,1);
xrect = [-1.5 1.5 0 1];
parity = 0;
elemtype = 0;
nodetype = 1;
[p,t] = squaremesh(ngrid,ngrid,parity,elemtype);
p(:,1) = xrect(1) + (xrect(2)-xrect(1))*p(:,1);
p(:,2) = xrect(3) + (xrect(4)-xrect(3))*p(:,2);
p(:,2) = loginc(p(:,2),2);
bndexpr = {'all(p(:,2)<min(p0(:,2))+1e-3)','all(p(:,1)>max(p0(:,1))-1e-3)', ...
           'all(p(:,2)>max(p0(:,2))-1e-3)','all(p(:,1)<min(p0(:,1))+1e-3)'};     
mesh = mkmesh(p,t,porder,bndexpr,elemtype,nodetype);

master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);

UDG = initu(mesh,{0;0;0});
UH = inituhat(master,mesh.elcon,UDG,app.ncu);

% HDG solver
[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);

% HDG postprocessing 
mesh1 = mkmesh_square(ngrid,ngrid,porder+1);
master1 = mkmaster(mesh1,2*(porder+1));
[master1,mesh1] = preprocess(master1,mesh1,hybrid);
UDGstar = postprocessnd(master,mesh,master1,mesh1,UDG);

figure(1); clf; scaplot(mesh,UDG(:,1,:),[],2); axis equal; axis tight;
%figure(2); clf; scaplot(mesh1,UDGstar(:,1,:),[],2); axis equal; axis tight;

app.adjoint = 1;
[VDG,VH] = hdg_solve(master,mesh,app,UDG,full(UH),[]);
figure(3);clf; scaplot(mesh,VDG(:,1,:),[],2); axis equal; axis tight;

