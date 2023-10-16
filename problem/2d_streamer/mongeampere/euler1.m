porder = 4;

nstage = 1;
torder = 1;
Mach   = 7;
aoa    = 0.0;
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
app.bcm  = [5,2,5,6];  
app.bcs  = [ui; ui; ui; ui];
app.bcd  = [1,1,1,1];  
app.bcv  = [0; 0; 0; 0];

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

app.fc_q = 1;
app.fc_u = 0;
app.tdep = false;
app.adjoint = 0;

mesh = mkmesh_square(31,31,porder,1,1,1,1,1);
mesh.p(:,1) = logdec(mesh.p(:,1),0.5);
mesh.dgnodes(:,1,:) = logdec(mesh.dgnodes(:,1,:),0.5);
mesh = mkmesh_halfcircle(mesh,1,3,4.7,pi/2,3*pi/2);

master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);
mesh = mkcgmesh(mesh);
mesh.ib = [];
mesh.in = 1:size(mesh.p2,1);
x = mesh.p2(:,1); y = mesh.p2(:,2);
mesh.ib = find(abs(x.^2 + y.^2)<1+1e-5);
mesh.in = setdiff( (1:size(mesh.p2,1))', mesh.ib);
[cgnodes,cgelcon,rowent2elem,colent2elem,cgent2dgent] = mkcgent2dgent(mesh.dgnodes(:,1:2,:),1e-8);
mesh.dist = tanh(meshdist(mesh,2)*20);

UDG = initu(mesh,{ui(1),ui(2),ui(3),ui(4),0,0,0,0,0,0,0,0});
UDG(:,2,:) = UDG(:,2,:).*tanh(meshdist(mesh,2)*4);
UDG(:,3,:) = UDG(:,3,:).*tanh(meshdist(mesh,2)*4);
TnearWall = pinf/(gam-1); % Tinf * (Twall/Tref-1) * exp(-10*dist) + Tinf;
UDG(:,4,:) = TnearWall + 0.5*(UDG(:,2,:).*UDG(:,2,:) + UDG(:,3,:).*UDG(:,3,:));
UH = inituhat(master,mesh.elcon,UDG,app.ncu);
SH = [];

mesh.dgnodes(:,3,:) = 0.05*mesh.dist;
[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,SH);

figure(1); clf; scaplot(mesh, eulereval(UDG,'M',1.4,7),[],2,1); 

[u,q,uhat,v,rho,iter] = hdgsolve_mah(master, mesh, [], UDG, 1e-3, mesh.dist);
[UDG1, UH1, ACG1, mine1, minf1, ming1, mesh1] = hdgsolve_avloop(master, mesh, app, q, UDG, UH, 0.05, 5);
figure(2); clf; scaplot(mesh1, eulereval(UDG1{5},'M',1.4,7),[],1,2); 

[u1,q1,uhat1,v1,rho1,iter1] = hdgsolve_mah(master, mesh, mesh1, UDG1{5}, 20e-3, []);
[elist, xi] = locatexinmesh(mesh1, q1, [], 1e-4);
UDG0 = evalfield(mesh1, UDG1{5}, elist, xi);
UDG0 = permute(reshape(UDG0, [master.npe mesh.ne 12]),[1 3 2]); 
UH0 = inituhat(master,mesh.elcon,UDG0,app.ncu);
[UDG2, UH2, ACG2, mine2, minf2, ming2, mesh2] = hdgsolve_avloop(master, mesh, app, q1, UDG0, UH0, 0.02, 5);
figure(3); clf; scaplot(mesh2, eulereval(UDG2{5},'M',1.4,7),[],1,2); axis on; axis equal; axis tight;

% mesht = mesh; mesht.dgnodes = q; figure(1); clf; meshplot(mesht,1); axis on; axis equal; axis tight;
% mesht = mesh; mesht.dgnodes = q1; figure(2); clf; meshplot(mesht,1); axis on; axis equal; axis tight;
% figure(3); clf; scaplot(mesh1, eulereval(UDG1{5},'M',1.4,7),[],1);  axis on; axis equal; axis tight;
% figure(4); clf; scaplot(mesht, eulereval(UDG0,'M',1.4,7),[],1);  axis on; axis equal; axis tight;

[u2,q2,uhat2,v2,rho2,iter2] = hdgsolve_mah(master, mesh, mesh2, UDG2{5}, 30e-3, [], 0.25);
[elist, xi] = locatexinmesh(mesh2, q2, [], 1e-4);
UDG0 = evalfield(mesh2, UDG2{5}, elist, xi);
UDG0 = permute(reshape(UDG0, [master.npe mesh.ne 12]),[1 3 2]); 
UH0 = inituhat(master,mesh.elcon,UDG0,app.ncu);
[UDG3, UH3, ACG3, mine3, minf3, ming3, mesh3] = hdgsolve_avloop(master, mesh, app, q2, UDG0, UH0, 0.02, 5);
figure(3); clf; scaplot(mesh3, eulereval(UDG3{5},'M',1.4,7),[],1,2); axis on; axis equal; axis tight;



