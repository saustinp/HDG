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

mesh = mkmesh_square(20,31,porder,1,1,1,1,1);
mesh.p(:,1) = logdec(mesh.p(:,1),1.5);
mesh.dgnodes(:,1,:) = logdec(mesh.dgnodes(:,1,:),1.5);
mesh = mkmesh_halfcircle(mesh,1,3,4.7,pi/2,3*pi/2);

master = mkmaster(mesh,2*porder);
mesh.dist=mshsize(mesh);
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
UH = inituhat(master,mesh.elcon,UDG,app.ncu);
SH = [];

app.fc_q = 1;
app.fc_u = 0;
app.tdep = false;
app.adjoint = 0;

mesh.dgnodes(:,3,:) = 0.05*mesh.dist;
[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,SH);

lambda01 = 0.04;
mesh.dgnodes(:,3,:) = lambda01*mesh.dist;
[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,SH);

S0 = 0.2;
kappa01 = 1.5;
eta = 0.8; m = 9;
lambda = ones(m,1)*lambda01;
for i = 2:m
    lambda(i) = lambda(i-1)*eta;
end
kappa = ones(m,1)*kappa01;
for i = 2:m
    kappa(i) = 1 + (kappa(i-1)-1)*eta;
end
[UDG1, UH1, ACG1, mine1, minf1, ming1] = avloop(master, mesh, app, UDG, UH, S0, lambda, kappa);


% for i = 1:length(UDG1)
%     figure(i); clf; scaplot(mesh, ACG1{i}(:,1,:),[],2,0); 
% end
% for i = 1:length(UDG1)
%     figure(i); clf; scaplot(mesh, UDG1{i}(:,1,:),[],2,0); 
% end
% for i = 1:length(UDG1)
%     figure(i); clf; scaplot(mesh, eulereval(UDG1{i},'M',1.4,7),[],2,0); 
% end

i = 7; nref=6; alpha = 100;  n = 40; m = 100; polydegree=6; 
div = divergence(UDG1{i}, 1.0);
divmax = max(div(:));
s = limiting(div,0,divmax/2,100,0);
sm = squeeze(mean(div,1));
idx = (sm > S0);
he = meshsize(mesh);
hs = he(:,idx);
hk =  (kappa(i)*min(hs(:)))^2;
s = cgpoisson(mesh, master, s, [hk 1.0]);        
s = s/max(s(:));
a = (s-S0).*(atan(alpha*(s-S0))/pi + 0.5) - atan(alpha)/pi + 0.5;    
S = reshape(a(mesh.t2'), master.npe, 1, mesh.ne);       
figure(1); clf; scaplot(mesh, S(:,1,:),[],2,0); 
[px1, sx1, idx1] = shockpoints(mesh, S, S0, S1, nref);
[pp1, pl1, pr1, pa1] = shockcells(px1, d, n, m);
x = squeeze(pa1(:,1,:))';
y = squeeze(pa1(:,2,:))';
[pn, un] = shocklocation(mesh, S, idx1, x, y, nref+2);
pn = pn(1:end-1,:);
pn = pn(2:end,:);

save result1.mat master mesh UDG UH UDG1 UH1 ACG1 mine1 minf1 ming1 lambda kappa pn un

i = 9;
figure(i); clf; scaplot(mesh2, eulereval(UDG2{i},'M',1.4,7),[],2,1); 
hold on;
plot(pn(:,1), pn(:,2), 'or');

% starting solution
% S0 = 0.25;
% kappa01 = 1.5;
% [UDG01, UH01, ACG01, mine01, minf01, ming01, lambda01] = avstart(master, mesh, app, UDG, UH, mesh.dgnodes(:,3,:), S0, kappa01);
% figure(1); clf; scaplot(mesh,ACG01,[],2,1); axis on; 
% figure(2); clf; scaplot(mesh,eulereval(UDG01,'M',1.4,7),[],2,1); axis on; 
 
% i = 8; 
% S = divergence(UDG1{i}, 1.0);
% [px1, sx1, idx1] = shockpoints(mesh, S, S0, S1, nref);
% [pp1, pl1, pr1, pa1] = shockcells(px1, d, n, m);
% x = squeeze(pa1(:,1,:))';
% y = squeeze(pa1(:,2,:))';
% [pn, un] = shocklocation(mesh, S, idx1, x, y, nref+2);

%save result1.mat master mesh UDG UH ACG UDG1 UH1 ACG1 mine1 minf1 ming1 lambda mine01 minf01 ming01 pn un

%  
% S1 = 0.4; nref = 6; d = 2; n = 40; m = 100; polydegree=6; 
% 
% i = 6; 
% S = divergence(UDG1{i}, 1.0);
% [px1, sx1, idx1] = shockpoints(mesh, S, S0, S1, nref);
% [pp1, pl1, pr1, pa1] = shockcells(px1, d, n, m);
% x = squeeze(pa1(:,1,:))';
% y = squeeze(pa1(:,2,:))';
% [pn1, un1] = shocklocation(mesh, S, idx1, x, y, nref+2);
% 
i = 8; 
S = divergence(UDG1{i}, 1.0);
[px1, sx1, idx1] = shockpoints(mesh, S, S0, S1, nref);
[pp1, pl1, pr1, pa1] = shockcells(px1, d, n, m);
x = squeeze(pa1(:,1,:))';
y = squeeze(pa1(:,2,:))';
[pn2, un2] = shocklocation(mesh, S, idx1, x, y, nref+2);
% 
% i = 10; 
% S = divergence(UDG1{i}, 1.0);
% [px1, sx1, idx1] = shockpoints(mesh, S, S0, S1, nref);
% [pp1, pl1, pr1, pa1] = shockcells(px1, d, n, m);
% x = squeeze(pa1(:,1,:))';
% y = squeeze(pa1(:,2,:))';
% [pn3, un3] = shocklocation(mesh, S, idx1, x, y, nref+2);
% 
% poly = polyfit(pn1(:,2) ,pn1(:,1), 6);
% y1 = linspace(pn1(1,2),pn1(end,2), 1000);
% x1 = polyval(poly,y1);
% 
% poly = polyfit(pn2(:,2) ,pn2(:,1), 6);
% y2 = linspace(pn2(1,2),pn2(end,2), 1000);
% x2 = polyval(poly,y2);
% 
% poly = polyfit(pn3(:,2) ,pn3(:,1), 6);
% y3 = linspace(pn3(1,2),pn3(end,2), 1000);
% x3 = polyval(poly,y3);
% 
% figure(1);clf; hold on;
% plot(x1,y1,'-r','LineWidth',2); 
% plot(x2,y2,'-b','LineWidth',2); 
% plot(x3,y3,'-g','LineWidth',2); 
% axis equal; axis tight;
% 

% figure(1); clf; hold on;
% scaplot(mesh,s,[0 0.5],2,1); %plot(px(:,1),px(:,2),'o');
% axis equal; axis tight;
% % 
% [pp, pl, pr, pa] = shockcells(px, d, n, m);
% x = squeeze(pa(:,1,:))';
% y = squeeze(pa(:,2,:))';
% [pn, un] = shocklocation(mesh, s, idx, x, y, nref+2);


% for i = 1:length(av)
%     figure(i); clf; scaplot(mesh, ACG1{i}(:,1,:),[],2,0); 
% end
% for i = 1:length(av)
%     figure(i); clf; scaplot(mesh, UDG1{i}(:,1,:),[],2,0); 
% end

% div = divergence(UDG, href);
% s = cgpoisson(mesh, master, div, [hk 1.0]);        
% a = (s-S0).*(atan(alpha*(s-S0))/pi + 0.5) - atan(alpha)/pi + 0.5;    
% a = reshape(a(mesh.t2'), master.npe, 1, mesh.ne);    
% figure(1); clf; scaplot(mesh, a,[],2,0); axis off;

% alpha = 100; beta = 0.5; href = 1.0; hk = 0.002;
% av = [1 1/2 1/4 1/8 1/16 1/32]*0.02;
% [UDG1, UH1, ACG1, mine1] = avloop(master, mesh, app, UDG, UH, av, href, hk, alpha, beta);

% alpha = 100;
% href = 1.0;
% div = divergence(UDG, href/5);
% s = ucg2udg(udg2ucg(div, cgent2dgent, rowent2elem), cgelcon);
% s = cgpoisson(mesh, master, s, [0.005 1.0]);        
% a = (s-0.1).*(atan(alpha*(s-0.1))/pi + 0.5) - atan(alpha)/pi + 0.5;
% a = reshape(a(mesh.t2'), master.npe, 1, mesh.ne);    
% figure(1); clf; scaplot(mesh, a,[],2,0); axis off;
% 
% alpha = 100; beta = 0.1; href = 0.2; hk = 0.002;
% av = [1 1/2 1/4 1/8 1/16 1/32]*0.1;
% [UDG1, UH1, ACG1, mine1] = avloop(master, mesh, app, UDG, UH, av, href, hk, alpha, beta);
% 
% for i = 1:length(av)
%     figure(i); clf; scaplot(mesh, ACG1{i}(:,1,:),[],2,0); 
% end
% for i = 1:length(av)
%     figure(i); clf; scaplot(mesh, UDG1{i}(:,1,:),[],2,0); 
% end
% 
% s = ucg2udg(udg2ucg(divergence(UDG1{5},href), cgent2dgent, rowent2elem), cgelcon);
% s0 = 0.1; s1 = 0.4; nref = 6; d = 2; n = 40; m = 100; polydegree=6; 
% [px, sx, idx] = shockpoints(mesh, s, s0, s1, nref);
% % figure(1); clf; hold on;
% % scaplot(mesh,s,[0 0.5],2,1); %plot(px(:,1),px(:,2),'o');
% % axis equal; axis tight;
% % 
% [pp, pl, pr, pa] = shockcells(px, d, n, m);
% x = squeeze(pa(:,1,:))';
% y = squeeze(pa(:,2,:))';
% [pn, un] = shocklocation(mesh, s, idx, x, y, nref+2);
% pn = pn(2:end-1,:);
% 
% figure(1); clf; hold on;
% scaplot(mesh,s,[0 0.5],2,1); plot(pn(:,1),pn(:,2),'o');
% axis equal; axis tight;
% 
% poly = polyfit(pn(:,2) ,pn(:,1), 6);
% y1 = linspace(pn(1,2),pn(end,2), 1000);
% x1 = polyval(poly,y1);
% figure(1);clf; hold on;
% scaplot(mesh,s,[0 0.5],2,1);
% plot(x1,y1,'-r','LineWidth',2); axis equal; axis tight;
% 
% figure(1); clf; hold on;
% scaplot(mesh,s,[0 0.5],2,1); plot(pp(:,1),pp(:,3),'o'); plot(pp(:,2),pp(:,4),'o');
% axis equal; axis tight;
% 
% figure(1); clf; hold on;
% plot(pp(:,1),pp(:,3),'o'); plot(pp(:,2),pp(:,4),'o');
% axis equal; axis tight;
% figure(1); clf; hold on;
% scaplot(mesh,s,[0 0.5],2,1);
% for i = 1:size(pp,1)
%     plot([pp(i,1) pp(i,2) pp(i,2) pp(i,1) pp(i,1)],[pp(i,3) pp(i,3) pp(i,4) pp(i,4) pp(i,3)],'-r');
%     plot(pa(:,1,i),pa(:,2,i),'o');
% end
% axis equal; axis tight;
% 
% figure(2); clf; scaplot(mesh,mesh.dgnodes(:,3,:),[],2,0); axis off; 
% axis equal; axis([-0.5 1.5 -1 1]);
% colorbar;

