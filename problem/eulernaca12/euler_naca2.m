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
app.source = 'source';
app.flux = 'flux';
app.fbou = 'fbou';
app.fhat = 'fhat';
app.hybrid = hybrid;
app.localsolve=0;
app.arg = {gam,epslm,tau};
app.bcm  = [2,1];  
app.bcs  = [ui; ui];
app.bcd  = [1,1];  % 2: Slip wall, 1: Far-field
app.bcv  = [0; 0];

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

SH = [];

app.fc_q = 1;
app.fc_u = 0;
app.tdep = false;
app.adjoint = 0;
[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);

figure(1); clf;
scaplot(mesh,UDG(:,1,:),[],2,0); axis off; 
axis equal; axis([-0.5 1.5 -1 1]);
colorbar;

