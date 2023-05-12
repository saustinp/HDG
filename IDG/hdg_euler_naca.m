
%NOTE: to keep this as close to Exasim as possible, I'm using Exasim preprocessing scripts when needed. 
addpath("exa_preprocessing_methods");
porder = 3;

nstage = 1;
torder = 1;
Mach   = 0.4;
aoa    = 1.5;
hybrid = 'hdg';

gam = 1.4;
epslm = 0.0;
Minf = Mach;                  % Infinity conditions
pinf = 1/(gam*Minf^2);
Re = inf;
Pr = 0.72;
alpha = aoa*pi/180;
tau = 2;

nd    = 2;
ntime = 30;
dt = 1e-3*2.^(0:ntime);
dt = repmat(dt,[nd 1]);
dt = dt(:);

ui = [ 1, cos(alpha), sin(alpha), 0.5+pinf/(gam-1)];

clear app;
app.source = 'euler_source';
app.flux = 'euler_flux';
app.fbou = 'euler_fbou';
app.fhat = 'euler_fhat';
app.hybrid = hybrid;
app.localsolve=0;
app.arg = {gam,epslm,tau};
app.bcm  = [2,1];  
app.bcs  = [ui; ui];
app.bcd  = [1,1];  % 2: Slip wall, 1: Far-field
app.bcv  = [0; 0];
app.ui = ui;

app.tdep = false;
app.wave = false;
app.alag = false;
app.flg_q = 1;
app.flg_p = 0;
app.flg_g = 0;

app.ndim = 2;
app.nch  = 2+app.ndim;                % Number of componets of UH
app.nc   = app.nch;                   % Number of componeents of UDG
app.ncu  = app.nch;                   % Number of components of U

app.time = [];
app.dtfc = [];
app.alpha = [];

mesh = mkmesh_naca12(porder);
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);

UDG = initu(mesh,{ui(1),ui(2),ui(3),ui(4)});
UH = inituhat(master,mesh.elcon,UDG,app.ncu);

SH = 0*UDG;

app.fc_q = 1;
app.fc_u = 0;
app.tdep = false;
app.adjoint = 0;
% [UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);
mesh.f2e = mkf2e(mesh.elcon, mesh.f, mesh.t2f, mesh.perm);
[nme, nbf] = exa_mkelemblocks(mesh.ne, 2000);
% I used a pared down version of Exasim's preprocessing script to 
% make sure that the faceblocks were being made correctly.
% Also, these scripts can construct the arrays (ent2ind1, rowe2f1, etc.)
% that are used for the element based face residual assembly. 
[app,mesh,~,dmd] = exa_preprocessing(app,mesh,master);
% 
nmf = app.fblks;
nn = [mesh.nd, app.nc, app.nc, mesh.nd, mesh.ne, nbf, mesh.ne, master.npv, master.ngv, mesh.nf, app.nbf(1), mesh.nf, master.npf, master.ngf];
% % cuda_check_grad(mesh, master, app, SH, UDG, UH, ui, app.arg, app.time, app.fc_q, app.fc_u, app.tdep, mesh.f2e, nn, nme, nmf)
cuda_check_inv(mesh, master, app, SH, UDG, UH, ui, app.arg, app.time, app.fc_q, app.fc_u, app.tdep, mesh.f2e, nn, nme, nmf)
% 
% figure(1); clf;
% scaplot(mesh,UDG(:,1,:),[],2,0); axis off; 
% axis equal; axis([-0.5 1.5 -1 1]);
% colorbar;

