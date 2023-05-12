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

mesh2 = naca0012shock(poly, poly2, porder);
master = mkmaster(mesh2,2*porder);
[master,mesh2] = preprocess(master,mesh2,hybrid);
mesh2 = mkcgmesh(mesh2);

dist = tanh(meshdist(mesh2,1)*30);
mesh2.dgnodes(:,3,:) = 0.01*dist;

UDG0 = initu(mesh2,{ui(1),ui(2),ui(3),ui(4),0,0,0,0,0,0,0,0});
UH0 = inituhat(master,mesh2.elcon,UDG0,app.ncu);

app.fc_q = 1;
app.fc_u = 0;
app.tdep = false;
app.adjoint = 0;

[UDG0,UH0] = hdg_solve(master,mesh2,app,UDG0,UH0,[]);
figure(2); clf;
scaplot(mesh2,UDG0(:,1,:),[],2,0); axis off; 
axis equal; axis([-0.5 1.5 -1 1]);
colorbar;

alpha = 100;
href = 0.1;
div = divergence(UDG0, href);
%s = ucg2udg(udg2ucg(div, cgent2dgent, rowent2elem), cgelcon);
s = cgpoisson(mesh2, master, div, [0.001 1.0]);        
a = (s-0.1).*(atan(alpha*(s-0.1))/pi + 0.5) - atan(alpha)/pi + 0.5;
a = reshape(a(mesh2.t2'), master.npe, 1, mesh2.ne);    
figure(2); clf; scaplot(mesh2, a,[],2,0); axis off;

mesh2.dgnodes(:,3,:) = max(0.0008*dist, 0.02.*a);
dist = tanh(meshdist(mesh2,1)*200); tanh(sqrt((mesh2.dgnodes(:,1,:)-1.0089).^2 + mesh2.dgnodes(:,2,:).^2)*120);
mesh2.dgnodes(:,3,:) = 0.0025.*a.*dist;
figure(2); clf; scaplot(mesh2, mesh2.dgnodes(:,3,:),[],2,0); axis off;
[UDG0,UH0] = hdg_solve(master,mesh2,app,UDG0,UH0,SH);
figure(1); clf; scaplot(mesh2, UDG0(:,1,:),[],2,0); axis off;


mesh2.dgnodes(:,3,:) = 0.0025.*a;
figure(1); clf; scaplot(mesh2, mesh2.dgnodes(:,3,:),[],2,0); axis off;
[UDG0,UH0] = hdg_solve(master,mesh2,app,UDG0,UH0,SH);
figure(1); clf; scaplot(mesh2, UDG0(:,1,:),[],2,0); axis off;

mesh2.dgnodes(:,3,:) = 0.0025.*a;

alpha = 100; beta = 0.1; href = 0.05; hk = 0.00005;
av = [1 1/2 1/4 1/8 1/16 1/32 1/64]*0.0025;
[UDG2, UH2, ACG2, mine2] = avloop(master, mesh2, app, UDG0, UH0, av, href, hk, alpha, beta);

x = linspace(-0.5, 1.5, 1200);
y = [-0.4 -0.2 -0.1 0.1 0.2 0.4 -0.15 0.15];
for i = 7:length(y)
    xy = [x(:) y(i)*ones(length(x),1)];
    rho1(:,i) = fieldatx(mesh,UDG1{3}(:,1,:),xy,7);    
    pres1(:,i) = fieldatx(mesh,eulereval(UDG1{3},'p',1.4,0.85),xy,7);    
    mach1(:,i) = fieldatx(mesh,eulereval(UDG1{3},'M',1.4,0.85),xy,7);    
end

%figure(4); clf; plot(x, mach1, '-');

for i = 7:length(y)
    xy = [x(:) y(i)*ones(length(x),1)];
    rho2(:,i) = fieldatx(mesh2,UDG2{6}(:,1,:),xy,7);    
    pres2(:,i) = fieldatx(mesh2,eulereval(UDG2{6},'p',1.4,0.85),xy,7);    
    mach2(:,i) = fieldatx(mesh2,eulereval(UDG2{6},'M',1.4,0.85),xy,7);    
end


%figure(1); clf; scaplot(mesh2, ACG2{end}(:,1,:),[],2,0); axis off;

% [master,mesh2] = preprocess(master,mesh2,hybrid);
% mesh2 = mkcgmesh(mesh2);
% mesh2.ib = [];
% mesh2.in = 1:size(mesh2.p2,1);
% [cgnodes,cgelcon,rowent2elem,colent2elem,cgent2dgent] = mkcgent2dgent(mesh2.dgnodes(:,1:2,:),1e-8);

% alpha = 100;
% href = 0.2;
% div = divergence(UDG0, href);
% s = ucg2udg(udg2ucg(div, cgent2dgent, rowent2elem), cgelcon);
% s = cgpoisson(mesh2, master, s, [0.001 1.0]);        
% a = (s-0.1).*(atan(alpha*(s-0.1))/pi + 0.5) - atan(alpha)/pi + 0.5;
% a = reshape(a(mesh2.t2'), master.npe, 1, mesh2.ne);    
% figure(2); clf; scaplot(mesh2, a,[],2,0); axis off;
% 
% %mesh2.dgnodes(:,3,:) = max(0.000*dist, 0.05.*a);
% mesh2.dgnodes(:,3,:) = 0.02.*a;
% figure(2); clf; scaplot(mesh2, mesh2.dgnodes(:,3,:),[],2,0); axis off;
% [UDG0,UH0] = hdg_solve(master,mesh2,app,UDG0,UH0,SH);
% figure(1); clf; scaplot(mesh2, UDG0(:,1,:),[],2,1); axis off;
% 
% alpha = 100; beta = 0.1; href = 0.2; hk = 0.001;
% av = [1 1/2 1/4 1/8 1/16 1/32]*0.02;
% [UDG01, UH01, ACG1, mine1] = avloop(master, mesh2, app, UDG0, UH0, av, href, hk, alpha, beta);
% 
% for i = 1:length(av)
%     figure(i); clf; scaplot(mesh2, ACG1{i}(:,1,:),[],2,0); axis([-0.2 1.2 -0.2 1.1]);
% end
% 
% s = ucg2udg(udg2ucg(divergence(UDG01{3},href), cgent2dgent, rowent2elem), cgelcon);
% s0 = 0.05; s1 = 0.35; nref = 6; d = 2; n = 22; m = 100; polydegree=6; 
% [px, sx, idx] = shockpoints(mesh2, s, s0, s1, nref);
% figure(1); clf; hold on;
% scaplot(mesh2,s,[0 2],2,1); plot(px(:,1),px(:,2),'o');
% axis equal; axis tight;
% 
% % upper shock
% ind = px(:,2) > 0.02;
% [pp, pl, pr, pa] = shockcells(px(ind,:), d, n, m);
% pp(1,3) = 0.046;
% figure(1); clf; hold on;
% scaplot(mesh2,s,[0 2],2,1); plot(pp(:,1),pp(:,3),'o'); plot(pp(:,2),pp(:,4),'o');
% axis equal; axis tight;
% 
% x = squeeze(pa(:,1,:))';
% y = squeeze(pa(:,2,:))';
% y(1,:) = 0.046;
% [pn, un] = shocklocation(mesh2, s, idx, x, y, nref);
% figure(1); clf; hold on;
% %scaplot(mesh2,s,[0 2],2,1); 
% plot(pp(:,1),pp(:,3),'o'); plot(pp(:,2),pp(:,4),'o'); plot(pn(:,1),pn(:,2),'*');
% axis equal; axis tight;
% 
% poly = polyfit(pn(:,2) ,pn(:,1), 3);
% y1 = linspace(pn(1,2),pn(end,2), 100);
% x1 = polyval(poly,y1);
% 
% figure(1); clf; hold on;
% mesh2plot(mesh2);
% for i = 1:size(pp,1)
%     plot([pp(i,1) pp(i,2) pp(i,2) pp(i,1) pp(i,1)],[pp(i,3) pp(i,3) pp(i,4) pp(i,4) pp(i,3)],'-b');
% end
% plot(pn(:,1),pn(:,2),'or');
% %plot(x1,y1,'-r');
% axis equal; axis tight;
% axis([-0.2 1.2 -0.2 0.8]);
% 
% figure(1); clf; hold on;
% mesh2plot(mesh2);
% plot(pn(:,1),pn(:,2),'or');
% plot(x1,y1,'-b','LineWidth',2);
% axis equal; axis tight;
% axis([-0.2 1.2 -0.2 0.8]);
% 
% figure(1); clf; hold on;
% scaplot(mesh2, UDG01{3}(:,1,:),[],2,1); axis on; colorbar off;
% plot(x1,y1,'-r','LineWidth',2);
% axis equal; axis tight;
% axis([-0.2 1.2 -0.2 0.8]);
% 
% % lower shock
% ind = px(:,2) < -0.02;
% [pp, pl, pr, pa] = shockcells(px(ind,:), d, n, m);
% pp(end,4) = -0.0605;
% figure(1); clf; hold on;
% scaplot(mesh2,s,[0 2],2,1); plot(pp(:,1),pp(:,3),'o'); plot(pp(:,2),pp(:,4),'o');
% axis equal; axis tight;
% 
% x = squeeze(pa(:,1,:))';
% y = squeeze(pa(:,2,:))';
% y(end,:) = -0.0605;
% [pn2, un2] = shocklocation(mesh2, s, idx, x, y, nref);
% figure(1); clf; hold on;
% plot(pp(:,1),pp(:,3),'o'); plot(pp(:,2),pp(:,4),'o'); plot(pn2(:,1),pn2(:,2),'*');
% axis equal; axis tight;
% 
% poly2 = polyfit(pn2(:,2) ,pn2(:,1), 3);
% y2 = linspace(pn2(1,2),pn2(end,2), 100);
% x2 = polyval(poly2,y2);
% 
% figure(1); clf; hold on;
% mesh2plot(mesh2);
% for i = 1:size(pp,1)
%     plot([pp(i,1) pp(i,2) pp(i,2) pp(i,1) pp(i,1)],[pp(i,3) pp(i,3) pp(i,4) pp(i,4) pp(i,3)],'-b');
% end
% plot(pn2(:,1),pn2(:,2),'or');
% plot(x2,y2,'-r');
% axis equal; axis tight;
% axis([-0.2 1.2 -0.2 0.8]);
% 
% 
% figure(1); clf; hold on;
% scaplot(mesh2, UDG01{3}(:,1,:),[],2,1); axis on; colorbar off;
% plot(x1,y1,'-r','LineWidth',2);
% plot(x2,y2,'-r','LineWidth',2);
% axis equal; axis tight;
% axis([-0.2 1.2 -0.2 0.8]);
