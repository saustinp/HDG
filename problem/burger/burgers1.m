porder = 2;
nstage = 1;
torder = 1;
ngrid = 12;
hybrid = 'hdg';

clear app;
app.source = 'source';
app.flux = 'flux';
app.fbou = 'fbou';
app.fhat = 'fhat';
app.hybrid = hybrid;
app.localsolve=1;
ntime  = 5;
dt = linspace(0.01,10,ntime);

kappa = 0.0;
c = [0.5,0.0];
tau = 1;
av = 0.0;

app.arg = {kappa,c,av,1/(ngrid*porder),tau};
app.bcm = [1;2;1;2];
app.bcs = [1;0;1;0];  
app.bcd  = [1,1,1,1]; 
app.bcv  = [0; 0; 0; 0];

app.tdep = false;
app.wave = false;
app.alag = false;
app.flg_q = 1;
app.flg_p = 0;
app.flg_g = 0;

app.nd = 2;
app.nch  = 1;                       % Number of componets of UH
app.nc   = app.nch*(app.nd+1);    % Number of componeents of UDG
app.ncu  = 1;

app.time = [];
app.dtfc = [];
app.alpha = [];
app.itmax = 1;
app.fc_q = 1;
app.fc_p = 0;
app.fc_u = 0;
app.tdep = false;
app.nc = 3;

% mesh and master
mkmesh_burger2;
mesh = mesh1;
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);

% initial solution
UDG = initu(mesh,{0.5;0;0});
UH=inituhat(master,mesh.elcon,UDG,1);

% CG discretization for AV field
[cgnodes,cgelcon,rowent2elem,colent2elem,cgent2dgent] = mkcgent2dgent(mesh.dgnodes(:,1:2,:),1e-8);
R = 0.1; wR = 0.1;
[nodeR,weightR] = nodelistcg(mesh,cgnodes,cgelcon,R,wR);

% HDG solver for constant viscosity field
mesh.dgnodes(:,3,:) = 0.0;
[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);
figure(1); clf;scaplot(mesh,UDG(:,1,:),[],2,1); axis equal; axis tight;

