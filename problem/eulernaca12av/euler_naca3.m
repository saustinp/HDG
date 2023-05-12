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

mesh = mkmesh_naca12(porder);
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

S0 = 0.2;
kappa01 = 2.5;
eta = 0.8; m = 10;
lambda = ones(m,1)*lambda01;
for i = 2:m
    lambda(i) = lambda(i-1)*eta;
end
kappa = ones(m,1)*kappa01;
for i = 2:m
    kappa(i) = 1 + (kappa(i-1)-1)*eta;
end
[UDG1, UH1, ACG1, mine1, minf1, ming1] = avloop2(master, mesh, app, UDG, UH, S0, lambda, kappa);

% for i = 1:length(UDG1)
%     figure(i); clf; scaplot(mesh, eulereval(UDG1{i},'p',1.4,0.8),[],2,0); 
% end
% 
% figure(2); clf; scaplot(mesh,divergence(UDG,1),[],2,0); axis off; 
% 
% S1=0.1;i = 12; nref=6; alpha = 100;  n = 40; m = 100; polydegree=3; d=2;
% div = divergence(UDG1{i}, 1.0);
% divmax = max(div(:));
% s = limiting(div,0,25,100,0);
% sm = squeeze(mean(div,1));
% idx = (sm > S1);
% he = meshsize(mesh);
% hs = he(:,idx);
% hk =  (kappa(i)*min(hs(:)))^2;
% s = cgpoisson(mesh, master, s, [hk 1.0]);        
% s = s/max(s(:));
% a = (s-S1).*(atan(alpha*(s-S1))/pi + 0.5) - atan(alpha)/pi + 0.5;    
% S = reshape(a(mesh.t2'), master.npe, 1, mesh.ne);       
% figure(1); clf; scaplot(mesh, S(:,1,:),[],2,0); 
% 
% % upper shock
% [px1, sx1, idx1] = shockpoints(mesh, S, S1, 2*S1, nref);
% ind = px1(:,2) > 0.02;
% [pp1, pl1, pr1, pa1] = shockcells(px1(ind,:), d, n, m);
% x = squeeze(pa1(:,1,:))';
% y = squeeze(pa1(:,2,:))';
% [pn1, un1] = shocklocation(mesh, S, idx1, x, y, nref);
% poly1 = polyfit(pn1(:,2) ,pn1(:,1), polydegree);
% y1 = linspace(pn1(1,2),pn1(end,2), 1000);
% x1 = polyval(poly1,y1);
% 
% figure(1); clf; scaplot(mesh, S(:,1,:),[],2,1); 
% hold on;
% plot(pn1(:,1),pn1(:,2),'o');
% plot(x1,y1,'-');

% lower shock
% n2 = 20; m2 = 100;
% [px2, sx2, idx2] = shockpoints(mesh, S, S1/2, S1, nref);
% ind = px2(:,2) < 0.02;
% [pp2, pl2, pr2, pa2] = shockcells(px2(ind,:), d, n2, m2);
% x = squeeze(pa2(:,1,:))';
% y = squeeze(pa2(:,2,:))';
% [pn2, un2] = shocklocation(mesh, S, idx2, x, y, nref);
% pn2 = pn2(3:end,:);
% pn2(end,:) = [];
% poly2 = polyfit(pn2(:,2) ,pn2(:,1), polydegree);
% y2 = linspace(pn2(1,2),pn2(end,2), 1000);
% x2 = polyval(poly2,y2);
% 
% figure(1); clf; scaplot(mesh, S(:,1,:),[],2,1); 
% hold on;
% plot(pn2(:,1),pn2(:,2),'o');
% plot(x2,y2,'-');

% save result1.mat master mesh UDG UH UDG1 UH1 ACG1 mine1 minf1 ming1 lambda kappa pn1 un1 poly1 pn2 un2 poly2
% 
% mesh2 = naca0012shock(poly1, poly2, porder);

% 
% [UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,SH);
% figure(2); clf;
% scaplot(mesh,UDG(:,1,:),[],2,0); axis off; 
% axis equal; axis([-0.5 1.5 -1 1]);
% colorbar;
% 
% [master,mesh] = preprocess(master,mesh,hybrid);
% mesh = mkcgmesh(mesh);
% mesh.ib = [];
% mesh.in = 1:size(mesh.p2,1);
% [cgnodes,cgelcon,rowent2elem,colent2elem,cgent2dgent] = mkcgent2dgent(mesh.dgnodes(:,1:2,:),1e-8);
% 
% alpha = 100;
% href = 0.2;
% div = divergence(UDG, href);
% s = ucg2udg(udg2ucg(div, cgent2dgent, rowent2elem), cgelcon);
% s = cgpoisson(mesh, master, s, [0.001 1.0]);        
% a = (s-0.1).*(atan(alpha*(s-0.1))/pi + 0.5) - atan(alpha)/pi + 0.5;
% a = reshape(a(mesh.t2'), master.npe, 1, mesh.ne);    
% figure(2); clf; scaplot(mesh, a,[],2,0); axis off;
% 
% %mesh.dgnodes(:,3,:) = max(0.000*dist, 0.05.*a);
% mesh.dgnodes(:,3,:) = 0.02.*a;
% figure(2); clf; scaplot(mesh, mesh.dgnodes(:,3,:),[],2,0); axis off;
% [UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,SH);
% figure(1); clf; scaplot(mesh, UDG(:,1,:),[],2,1); axis off;
% 
% alpha = 100; beta = 0.1; href = 0.2; hk = 0.001;
% av = [1 1/2 1/4 1/8 1/16 1/32]*0.02;
% [UDG1, UH1, ACG1, mine1] = avloop(master, mesh, app, UDG, UH, av, href, hk, alpha, beta);
% 
% for i = 1:length(av)
%     figure(i); clf; scaplot(mesh, ACG1{i}(:,1,:),[],2,0); axis([-0.2 1.2 -0.2 1.1]);
% end
% 
% s = ucg2udg(udg2ucg(divergence(UDG1{3},href), cgent2dgent, rowent2elem), cgelcon);
% s0 = 0.05; s1 = 0.35; nref = 6; d = 2; n = 22; m = 100; polydegree=6; 
% [px, sx, idx] = shockpoints(mesh, s, s0, s1, nref);
% figure(1); clf; hold on;
% scaplot(mesh,s,[0 2],2,1); plot(px(:,1),px(:,2),'o');
% axis equal; axis tight;
% 
% % upper shock
% ind = px(:,2) > 0.02;
% [pp, pl, pr, pa] = shockcells(px(ind,:), d, n, m);
% pp(1,3) = 0.046;
% figure(1); clf; hold on;
% scaplot(mesh,s,[0 2],2,1); plot(pp(:,1),pp(:,3),'o'); plot(pp(:,2),pp(:,4),'o');
% axis equal; axis tight;
% 
% x = squeeze(pa(:,1,:))';
% y = squeeze(pa(:,2,:))';
% y(1,:) = 0.046;
% [pn, un] = shocklocation(mesh, s, idx, x, y, nref);
% figure(1); clf; hold on;
% %scaplot(mesh,s,[0 2],2,1); 
% plot(pp(:,1),pp(:,3),'o'); plot(pp(:,2),pp(:,4),'o'); plot(pn(:,1),pn(:,2),'*');
% axis equal; axis tight;
% 
% poly = polyfit(pn(:,2) ,pn(:,1), 3);
% y1 = linspace(pn(1,2),pn(end,2), 100);
% x1 = polyval(poly,y1);
% 
% figure(1); clf; hold on;
% meshplot(mesh);
% for i = 1:size(pp,1)
%     plot([pp(i,1) pp(i,2) pp(i,2) pp(i,1) pp(i,1)],[pp(i,3) pp(i,3) pp(i,4) pp(i,4) pp(i,3)],'-b');
% end
% plot(pn(:,1),pn(:,2),'or');
% %plot(x1,y1,'-r');
% axis equal; axis tight;
% axis([-0.2 1.2 -0.2 0.8]);
% 
% figure(1); clf; hold on;
% meshplot(mesh);
% plot(pn(:,1),pn(:,2),'or');
% plot(x1,y1,'-b','LineWidth',2);
% axis equal; axis tight;
% axis([-0.2 1.2 -0.2 0.8]);
% 
% figure(1); clf; hold on;
% scaplot(mesh, UDG1{3}(:,1,:),[],2,1); axis on; colorbar off;
% plot(x1,y1,'-r','LineWidth',2);
% axis equal; axis tight;
% axis([-0.2 1.2 -0.2 0.8]);
% 
% % lower shock
% ind = px(:,2) < -0.02;
% [pp, pl, pr, pa] = shockcells(px(ind,:), d, n, m);
% pp(end,4) = -0.0605;
% figure(1); clf; hold on;
% scaplot(mesh,s,[0 2],2,1); plot(pp(:,1),pp(:,3),'o'); plot(pp(:,2),pp(:,4),'o');
% axis equal; axis tight;
% 
% x = squeeze(pa(:,1,:))';
% y = squeeze(pa(:,2,:))';
% y(end,:) = -0.0605;
% [pn2, un2] = shocklocation(mesh, s, idx, x, y, nref);
% figure(1); clf; hold on;
% plot(pp(:,1),pp(:,3),'o'); plot(pp(:,2),pp(:,4),'o'); plot(pn2(:,1),pn2(:,2),'*');
% axis equal; axis tight;
% 
% poly2 = polyfit(pn2(:,2) ,pn2(:,1), 3);
% y2 = linspace(pn2(1,2),pn2(end,2), 100);
% x2 = polyval(poly2,y2);
% 
% figure(1); clf; hold on;
% meshplot(mesh);
% for i = 1:size(pp,1)
%     plot([pp(i,1) pp(i,2) pp(i,2) pp(i,1) pp(i,1)],[pp(i,3) pp(i,3) pp(i,4) pp(i,4) pp(i,3)],'-b');
% end
% plot(pn(:,1),pn(:,2),'or');
% plot(x2,y2,'-r');
% axis equal; axis tight;
% axis([-0.2 1.2 -0.2 0.8]);
% 
% 
% figure(1); clf; hold on;
% scaplot(mesh, UDG1{3}(:,1,:),[],2,1); axis on; colorbar off;
% plot(x1,y1,'-r','LineWidth',2);
% plot(x2,y2,'-r','LineWidth',2);
% axis equal; axis tight;
% axis([-0.2 1.2 -0.2 0.8]);
