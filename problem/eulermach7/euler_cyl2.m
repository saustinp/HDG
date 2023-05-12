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
[cgnodes,cgelcon,rowent2elem,colent2elem,cgent2dgent] = mkcgent2dgent(mesh2.dgnodes(:,1:2,:),1e-8);
dist = tanh(meshdist(mesh2,2)*40);

UDG02 = initu(mesh2,{ui(1),ui(2),ui(3),ui(4),0,0,0,0,0,0,0,0});
UH02 = inituhat(master,mesh2.elcon,UDG02,app.ncu);
SH = [];

app.fc_q = 1;
app.fc_u = 0;
app.tdep = false;
app.adjoint = 0;

mesh2.dgnodes(:,3,:) = 0.05*dist;
[UDG02,UH02] = hdg_solve(master,mesh2,app,UDG02,UH02,SH);

mesh2.dgnodes(:,3,:) = 0.025*dist;
[UDG02,UH02] = hdg_solve(master,mesh2,app,UDG02,UH02,SH);

% figure(2); clf; scaplot(mesh2,UDG02(:,1,:),[],2,0); axis off; 
% hold on; 
% plot(pn(:,1),pn(:,2),'-r','LineWidth',3);

mesh2.dgnodes(:,3,:) = 0.0125*dist;
[UDG02,UH02] = hdg_solve(master,mesh2,app,UDG02,UH02,SH);
figure(2); clf; scaplot(mesh2,UDG02(:,1,:),[],2,0); axis off; 

alpha = 100; beta = 0.1; href = 0.1; hk = 0.001;
av = [1 1/2 1/4 1/8 1/16 1/32 1/64 1/128]*0.1;
[UDG2, UH2, ACG2, mine2] = avloop(master, mesh2, app, UDG02, UH02, av, href, hk, alpha, beta);

for i = 1:length(av)
    figure(i); clf; scaplot(mesh2, ACG2{i}(:,1,:),[],2,0); 
end
for i = 1:length(av)
    figure(i); clf; scaplot(mesh2, UDG2{i}(:,1,:),[],2,0); 
end

% mesh2 = mkmesh_square(20,51,porder,1,1,1,1,1);
% p = mesh2.p;
% p(:,1) = loginc(logdec(p(:,1),3),7.5);
% %mesh2.dgnodes(:,1,:) = logdec(mesh2.dgnodes(:,1,:),1.5);
% 
% a=1; % inner radius
% t1 = pi/2;
% t2 = 3*pi/2;
% t = t2 + p(:,2)*(t1-t2);
% 
% rn = sqrt(pn(:,1).^2 + pn(:,2).^2);
% tn = atan(pn(:,2)./pn(:,1)) + pi;
% % nx = cos(tn);
% % ny = sin(tn);
% m = 51;
% tm = t2 + linspace(0,1,m)'*(t1-t2);
% nx = cos(tm);
% ny = sin(tm);
% 
% polyn = polyfit(tn, rn, 6);
% 
% dm = polyval(polyn, tm); 
% xm = dm.*cos(tm);
% ym = dm.*sin(tm);
% ds = [0 0.001 0.003 0.007 0.016 0.032 0.064 0.128];
% ds = [-ds(end:-1:2) ds];
% ds = [ds 0.21 0.4 0.7 1.1734];
% xn = zeros(m, length(ds));
% yn = 0*xn;
% for i = 1:length(ds)
%     xn(:,i) = xm + ds(i)*nx;
%     yn(:,i) = ym + ds(i)*ny;
% end
% 
% polym = polyfit(tm, sqrt(xn(:,1).^2 + yn(:,1).^2), 6);
% d = polyval(polym, t); % outer radius
% r = d + p(:,1).*(a-d);
% pnew = 0*p;
% pnew(:,1) = r.*cos(t);
% pnew(:,2) = r.*sin(t);
% 
% figure(3); clf; hold on;
% simpplot(pnew, mesh2.t); 
% plot(xm,ym,'-or');
% for i = 1:length(ds)
%     plot(xn(:,i),yn(:,i),'-','LineWidth',3);
% end
% axis equal; axis tight;
% 
% 
% 
% % %dn = sqrt(xn(:,1).^2 + yn(:,1).^2);
% % %rn = dn + p(:,1).*(a-dn);
% % 
% % d = polyval(polyn, t); % outer radius
% % r = d + p(:,1).*(a-d);
% % pnew = 0*p;
% % pnew(:,1) = r.*cos(t);
% % pnew(:,2) = r.*sin(t);
% % 
% % poly = polyfit(pn(:,2) ,pn(:,1), 6);
% % y1 = linspace(pn(1,2),pn(end,2), 1000);
% % x1 = polyval(poly,y1);
% % 
% % figure(1); clf; hold on;
% % simpplot(pnew, mesh2.t); 
% % %simpplot(mesh.p, mesh.t); 
% % %plot(x1,y1,'o');
% % axis on;
% % 
% % figure(2); clf; hold on;
% % simpplot(mesh.p, mesh.t); 
% % %plot(x1,y1,'o');
% % axis on;
% 
% 
