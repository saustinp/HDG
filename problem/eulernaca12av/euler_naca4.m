porder = 4;

nstage = 1;
torder = 1;
Mach   = 0.8;
aoa    = 1.5;
hybrid = 'hdg';

gam = 1.4;
epslm = 0.0;
Minf = Mach;                  % Infinity conditions
pinf = 1/(gam*Minf^2);
Re = inf;
Pr = 0.72;
alpha = aoa*pi/180;
tau = 1;

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
app.localsolve = 1;
app.arg = {gam,Minf,epslm,tau};
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
app.nc   = app.nch*3;                   % Number of componeents of UDG
app.ncu  = app.nch;                   % Number of components of U

app.time = [];
app.dtfc = [];
app.alpha = [];

%mesh = mkmesh_naca12(porder);
mesh = naca0012shock3(polya, polyb, porder);
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);
mesh = mkcgmesh(mesh);
x = mesh.p2(:,1); y = naca12(x);
ib1 = find(abs(mesh.p2(:,2)-y)<1e-6);
ib2 = find(abs(mesh.p2(:,2)+y)<1e-6);
mesh.ib = unique([ib1; ib2]);
mesh.in = setdiff( (1:size(mesh.p2,1))', mesh.ib);
mesh.dist = tanh(meshdist(mesh,1)*30);

UDG = initu(mesh,{ui(1),ui(2),ui(3),ui(4),0,0,0,0,0,0,0,0});
UH = inituhat(master,mesh.elcon,UDG,app.ncu);
SH = [];

app.fc_q = 1;
app.fc_u = 0;
app.tdep = false;
app.adjoint = 0;

lambda01 = 0.01;
mesh.dgnodes(:,3,:) = lambda01*mesh.dist;
[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,SH);

lambda01 = 0.005;
mesh.dgnodes(:,3,:) = lambda01*mesh.dist;
[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,SH);

lambda01 = 0.0025;
mesh.dgnodes(:,3,:) = lambda01*mesh.dist;
[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,SH);

[UDG1, UH1, ACG1, mine1, minf1, ming1] = avloop2(master, mesh, app, UDG, UH, 0.1/4, lambda01, 5);
[UDG1, UH1, ACG1, mine1, minf1, ming1] = avloop2(master, mesh, app, UDG1{1}, UH1{1}, 0.1/4, lambda01, 5);
UDG02=UDG1{1};
UH02=UH1{1};
 
S0 = 0.1;
kappa01 = 5;
eta = 0.8; m = 14;
lambda = ones(m,1)*lambda01;
for i = 2:m
    lambda(i) = lambda(i-1)*eta;
end
kappa = ones(m,1)*kappa01;
for i = 2:m
    kappa(i) = 1 + (kappa(i-1)-1)*eta;
end
[UDG1, UH1, ACG1, mine1, minf1, ming1] = avloop2(master, mesh, app, UDG02, UH02, S0, lambda, kappa);


% figure(1); clf; scaplot(mesh, eulereval(UDG1{14},'p',1.4,0.8),[],2,0); 
% porder = 4;
% 
% nstage = 1;
% torder = 1;
% Mach   = 0.8;
% aoa    = 1.5;
% hybrid = 'hdg';
% 
% gam = 1.4;
% epslm = 0.0;
% Minf = Mach;                  % Infinity conditions
% pinf = 1/(gam*Minf^2);
% Re = inf;
% Pr = 0.72;
% alpha = aoa*pi/180;
% tau = 1;
% 
% nd    = 2;
% ntime = 30;
% dt = 1e-3*2.^(0:ntime);
% dt = repmat(dt,[nd 1]);
% dt = dt(:);
% 
% ui = [ 1, cos(alpha), sin(alpha), 0.5+pinf/(gam-1)];
% 
% clear app;
% app.source = 'source';
% app.flux = 'flux';
% app.fbou = 'fbou';
% app.fhat = 'fhat';
% app.hybrid = hybrid;
% app.localsolve = 1;
% app.arg = {gam,Minf,epslm,tau};
% app.bcm  = [2,1];  
% app.bcs  = [ui; ui];
% app.bcd  = [1,1];  % 2: Slip wall, 1: Far-field
% app.bcv  = [0; 0];
% 
% app.tdep = false;
% app.wave = false;
% app.alag = false;
% app.flg_q = 1;
% app.flg_p = 0;
% app.flg_g = 0;
% 
% app.ndim = 2;
% app.nch  = 2+app.ndim;                % Number of componets of UH
% app.nc   = app.nch*3;                   % Number of componeents of UDG
% app.ncu  = app.nch;                   % Number of components of U
% 
% app.time = [];
% app.dtfc = [];
% app.alpha = [];
% 
% %mesh2 = naca0012shock2(poly1, poly2, porder);
% mesh2 = naca0012shock(polya, polyb, porder);
% master = mkmaster(mesh2,2*porder);
% [master,mesh2] = preprocess(master,mesh2,hybrid);
% mesh2 = mkcgmesh(mesh2);
% x = mesh2.p2(:,1); y = naca12(x);
% ib1 = find(abs(mesh2.p2(:,2)-y)<1e-5);
% ib2 = find(abs(mesh2.p2(:,2)+y)<1e-5);
% mesh2.ib = unique([ib1; ib2]);
% mesh2.in = setdiff( (1:size(mesh2.p2,1))', mesh2.ib);
% mesh2.dist = tanh(meshdist(mesh2,1)*30);
% 
% figure(1);clf;
% meshplot(mesh2,[0 1]); hold on;
% plot(mesh2.p2(mesh2.ib,1),mesh2.p2(mesh2.ib,2),'*'); axis equal; axis tight
% 
% UDG02 = initu(mesh2,{ui(1),ui(2),ui(3),ui(4),0,0,0,0,0,0,0,0});
% UH02 = inituhat(master,mesh2.elcon,UDG02,app.ncu);
% SH = [];
% 
% app.fc_q = 1;
% app.fc_u = 0;
% app.tdep = false;
% app.adjoint = 0;
% 
% lambda02 = 0.01;
% mesh2.dgnodes(:,3,:) = lambda02*mesh2.dist;
% [UDG02,UH02] = hdg_solve(master,mesh2,app,UDG02,UH02,SH);
% 
% lambda02 = 0.005;
% mesh2.dgnodes(:,3,:) = lambda02*mesh2.dist;
% [UDG02,UH02] = hdg_solve(master,mesh2,app,UDG02,UH02,SH);
% 
% lambda02 = 0.0025;
% mesh2.dgnodes(:,3,:) = lambda02*mesh2.dist;
% [UDG02,UH02] = hdg_solve(master,mesh2,app,UDG02,UH02,SH);
% 
% figure(1); clf; scaplot(mesh2, eulereval(UDG02,'p',1.4,0.8),[],2,0); 
% 
% S0 = 0.2;
% kappa02 = 10;
% eta = 0.7; m = 16;
% lambda2 = ones(m,1)*lambda02;
% for i = 2:m
%     lambda2(i) = lambda2(i-1)*eta;
% end
% kappa2 = ones(m,1)*kappa02;
% for i = 2:m
%     kappa2(i) = 1 + (kappa2(i-1)-1)*eta;
% end
% [UDG2, UH2, ACG2, mine2, minf2, ming2] = avloop(master, mesh2, app, UDG02, UH02, S0/10, lambda2, kappa2);
% 
% 
% for i = 1:length(UDG2)
%     figure(i); clf; scaplot(mesh2, eulereval(UDG2{i},'p',1.4,0.8),[],2,0); 
% end
% 
% 
% S0 = 0.2;
% kappa02 = 5;
% eta = 0.8; m = 16;
% lambda2 = ones(m,1)*lambda02/10;
% for i = 2:m
%     lambda2(i) = lambda2(i-1)*eta;
% end
% kappa2 = ones(m,1)*kappa02;
% for i = 2:m
%     kappa2(i) = 1 + (kappa2(i-1)-1)*eta;
% end
% [UDG2, UH2, ACG2, mine2, minf2, ming2] = avloop(master, mesh2, app, UDG01, UH01, S0, lambda2, kappa2);
% 
