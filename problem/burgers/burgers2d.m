setapplicationpath('FM/viscburgers');

porder = 6;
h  = 1/12;
nstage = 1;
torder = 1;
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

kappa = 0;
c = [0.5,0.5];
tau = 1;

app.arg = {kappa,c,tau};
app.bcm = [1;1;1;1];
app.bcs = [0;0;0;0]; %[1,0,0;1,0,0];
app.bcd  = [1,1,1,1];  % 2: Slip wall, 1: Far-field
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

mesh = mkmesh_unstructuredsquare(h,porder);
master = mkmaster(mesh,2*porder);

[master,mesh] = preprocess(master,mesh,hybrid);

UDG = initu(mesh,{0;0;0});
UH=inituhat(master,mesh.elcon,UDG,1);


% HDG solver
app.fc_q = 1;
app.fc_p = 0;
app.fc_u = 0;
app.tdep = false;
app.nc = 3;
[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);
figure(2); clf;scaplot(mesh,UDG(:,1,:),[],2,1); axis equal; axis tight;

% % HDG postprocessing 
% mesh1 = mkmesh_square(ngrid,ngrid,porder+1);
% master1 = mkmaster(mesh1,2*(porder+1));
% [master1,mesh1] = preprocess(master1,mesh1,hybrid);
% UDGstar = postprocessnd(master,mesh,master1,mesh1,UDG);
% 
% figure(1); clf;scaplot(mesh,UDG(:,1,:),[],2,1); axis equal; axis tight;
% figure(2); clf;scaplot(mesh1,UDGstar(:,1,:),[],2,1); axis equal; axis tight;
