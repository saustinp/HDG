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
app.bcm  = [5,2,6];  
app.bcs  = [ui; ui; ui];
app.bcd  = [1,1,1];  
app.bcv  = [0; 0; 0];

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

mesh2 = mkmesh_cylshock2(porder, pn);
master = mkmaster(mesh2,2*porder);
mesh2.dist=mshsize(mesh2);
[master,mesh2] = preprocess(master,mesh2,hybrid);
mesh2 = mkcgmesh(mesh2);
mesh2.ib = [];
mesh2.in = 1:size(mesh2.p2,1);
x = mesh2.p2(:,1); y = mesh2.p2(:,2);
mesh2.ib = find(abs(x.^2 + y.^2)<1+1e-5);
mesh2.in = setdiff( (1:size(mesh2.p2,1))', mesh2.ib);
[cgnodes,cgelcon,rowent2elem,colent2elem,cgent2dgent] = mkcgent2dgent(mesh2.dgnodes(:,1:2,:),1e-8);
dist = tanh(meshdist(mesh2,2)*20);

UDG02 = initu(mesh2,{ui(1),ui(2),ui(3),ui(4),0,0,0,0,0,0,0,0});
UH02 = inituhat(master,mesh2.elcon,UDG02,app.ncu);
SH = [];

app.fc_q = 1;
app.fc_u = 0;
app.tdep = false;
app.adjoint = 0;

mesh2.dgnodes(:,3,:) = 0.05*dist;
[UDG02,UH02] = hdg_solve(master,mesh2,app,UDG02,UH02,SH);

mesh2.dgnodes(:,3,:) = 0.02*dist;
[UDG02,UH02] = hdg_solve(master,mesh2,app,UDG02,UH02,SH);

% mesh2.dgnodes(:,3,:) = 0.0125*dist;
% [UDG02,UH02] = hdg_solve(master,mesh2,app,UDG02,UH02,SH);
% figure(2); clf; scaplot(mesh2,eulereval(UDG02,'M',1.4,7),[],2,1); axis on; 

S0 = 0.2; lambda02 = 0.02; 
kappa02 = 5;
eta = 0.8; m = 16;
lambda2 = ones(m,1)*lambda02;
for i = 2:m
    lambda2(i) = lambda2(i-1)*eta;
end
kappa2 = ones(m,1)*kappa02;
for i = 2:m
    kappa2(i) = 1 + (kappa2(i-1)-1)*eta;
end
[UDG2, UH2, ACG2, mine2, minf2, ming2] = avloop(master, mesh2, app, UDG02, UH02, S0, lambda2, kappa2);

%save result2.mat mesh2 UDG02 UH02 UDG2 UH2 ACG2 mine2 minf2 ming2 lambda2 kappa2


% %lambda = 0.003;
% eta = 0.8;
% av = ones(8,1)*0.0125;
% for i = 2:length(av)
%     av(i) = av(i-1)*eta;
% end
% [UDG2, UH2, ACG2, mine2, minf2, ming2] = avloop(master, mesh2, app, UDG01, UH01, av, S0, Chk);

% starting solution
% S0 = 0.5; Chk = 2;
% [UDG02, UH02, ACG02, mine02, minf02, ming02, lambda2] = avstart(master, mesh2, app, UDG01, UH01, mesh2.dgnodes(:,3,:), S0, Chk);
% figure(1); clf; scaplot(mesh2,eulereval(UDG01,'M',1.4,7),[],2,1); axis on; 
% figure(2); clf; scaplot(mesh2,eulereval(UDG02,'M',1.4,7),[],2,1); axis on; 
% figure(3); clf; scaplot(mesh2,ACG02,[],2,1); axis on; 

% for i = 1:length(UDG2)
%     figure(i); clf; scaplot(mesh2, UDG2{i}(:,1,:),[],2,0); 
% end
for i = 1:length(UDG2)
    figure(i); clf; scaplot(mesh2, eulereval(UDG2{i},'M',1.4,7),[],2,0); 
end

% figure(1); clf; scaplot(mesh2,eulereval(UDG2{2},'M',1.4,7),[],2,1); axis off; 

% [UDG3, UH3, ACG3, mine3, minf3, ming3] = avloop(master, mesh2, app, UDG02, UH02, 0.003, S0);
% figure(1); clf; scaplot(mesh2,eulereval(UDG3{1},'M',1.4,7),[],2,1); axis on; 
% figure(2); clf; scaplot(mesh2, ACG3{1}(:,1,:),[],2,0); axis on;
% figure(2); clf; scaplot(mesh2, divergence(UDG02, 1.0),[],2,0); 


% i = 5; 
% S = divergence(UDG2{i}, 1.0);
% figure(i); clf; scaplot(mesh2, S,[],2,0); 
% [px1, sx1, idx1] = shockpoints(mesh2, S, S0, S1, nref);
% [pp1, pl1, pr1, pa1] = shockcells(px1, d, n, m);
% x = squeeze(pa1(:,1,:))';
% y = squeeze(pa1(:,2,:))';
% [pn4, un4] = shocklocation(mesh2, S, idx1, x, y, nref+2);

%figure(1); clf; plot(x(:), y(:),'o');

% alpha = 100; beta = 0.1; href = 0.1; hk = 0.001;
% av = [1 1/2 1/4 1/8 1/16 1/32 1/64 1/128]*0.1;
% [UDG2, UH2, ACG2, mine2] = avloop(master, mesh2, app, UDG02, UH02, av, href, hk, alpha, beta);
% 
% for i = 1:length(av)
%     figure(i); clf; scaplot(mesh2, ACG2{i}(:,1,:),[],2,0); 
% end
% for i = 1:length(av)
%     figure(i); clf; scaplot(mesh2, UDG2{i}(:,1,:),[],2,0); 
% end

