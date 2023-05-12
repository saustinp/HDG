setapplicationpath('FM/ns');

m      = 21;
n      = 31;
porder = 4;
nstage = 1;
torder = 1;
hybrid = 'hdg';

gam = 1.4;
epslm = 0.1;
Minf = 0.2;                  % Infinity conditions
pinf = 1/(gam*Minf^2);
alpha = 0*pi/180;
Re = 100;
Pr = 0.72;
tau = 3.0;

ntime = 4;
%dt = linspace(0.01,1,ntime);
dt = [1e-4 1e-3 1e-2 1e-1]*1e-1;

ui = [ 1, cos(alpha), sin(alpha), 0.5+pinf/(gam-1)];
%uip = [pinf, pinf, pinf, pinf];

clear app;
app.source = 'source';
app.flux = 'flux';
app.fbou = 'fbou';
app.fhat = 'fhat';

app.getdqdg = 1;
app.denseblock = 0;
app.hybrid = hybrid;
app.localsolve = 1;
app.arg = {gam,epslm,Re,Pr,Minf,tau};
app.bcm  = [2,1];  % 2: Wall, 1: Far-field
app.bcs  = [ui;ui];

app.bcd  = [1,3];  % 2: Wall, 1: Far-field
app.bcv  = [ui;ui];

app.wave = false;
app.tdep = true;
app.alag = false;
app.flg_q = 1;
app.flg_p = 0;
app.flg_g = 0;
app.fc_q = 1;
app.fc_u = 1;
app.fc_p = 0;

app.np   = 2;
app.nd   = 2;
app.ncu  = 2+app.nd;               % Number of components of U
app.nch  = app.ncu;                % Number of componets of UH
app.nq   = app.ncu*app.nd;         % Number of componets of Q
app.nc   = app.ncu+app.nq;         % Number of componeents of UDG

app.time = [];
app.dtfc = [];
app.alpha = [];

mesh = mkmesh_trefftz(m,n,porder,[0.05,0.05,1.98]);
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);

UDG = initu(mesh,{ui(1),ui(2),ui(3),ui(4); 0,0,0,0; 0,0,0,0});
UH = inituhat(master,mesh.elcon,UDG,app.ncu);
% UDG = UDG+0.2*rand(size(UDG));
% UH = UH+0.2*rand(size(UH));

% time = 0;
% figure(1); clf;
% for itime = 1:ntime
%     fprintf('Timestep :  %d\n', itime);
%     
%     [UDG,UH] = hdg_solve_dirk(master,mesh,app,UDG,UH,[],time,dt(itime),nstage,torder);
%     time = time + dt(itime); 
%     
%     scaplot(mesh,eulereval(UDG(:,1:app.nch,:),'M',gam),[],1,1); axis off; %axis equal; axis tight;
%     tm=sprintf('t* = %.3f',time);
%     title(tm,'Color','white','FontSize',16,'FontName','Courier');
% end
% 

fprintf('\nSteady State Solve\n');

app.fc_q = 1;
app.fc_u = 0;
app.tdep = false;
app.adjoint = 0;
[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);

figure(1); clf; scaplot(mesh,eulereval(UDG(:,1:app.nch,:),'M',gam),[],2); axis off; %axis equal; axis tight;
% %     tm=sprintf('t* = %.3f',time);
% %     title(tm,'Color','white','FontSize',16,'FontName','Courier');

app.adjoint = 1;
[VDG,VH] = hdg_solve(master,mesh,app,UDG,UH,[]);
figure(2); clf; scaplot(mesh,VDG(:,1,:),[],2); axis off;

% figure(1);clf;
% scaplot(mesh,eulereval(UDG(:,1:app.nch,:),'M',gam),[],2,0); 
% axis off;
% axis([-3 4 -2.5 2.5]);
% colorbar('FontSize',14);
% 

% load hdgp2.mat;
% [Cp,Cf,x]=getsurfacedata(master,mesh,app,UDG,UH,2);
% load edgp2.mat;
% [Cp2,Cf2,x2]=getsurfacedata(master,mesh,app,UDG,UH,2);
% 
% figure(1); clf;
% set(axes, 'FontSize', 18, 'LineWidth', 1.0, 'TickLength', [0.015 0.015], 'GridLineStyle', '-', ...
%      'XTick', [0 0.2 0.4 0.6 0.8 1],'YTick', [-1.0:0.2:0.4], 'YDir', 'Reverse');
% plot(x,Cp,'-b',x2,Cp2,'--r','LineWidth',1.5); 
% xlabel('$x/c$','FontSize', 24, 'Interpreter', 'latex');
% ylabel('$C_p$','FontSize', 24, 'Interpreter', 'latex');
% legend({'HDG', 'EDG'},'Location','SE','Interpreter','latex');
% box on;
% grid on;
% axis on;
% axis([0 1 -1 0.4]);
% print -depsc nstrefftzcpp2.eps
% 
% figure(2); clf;
% set(axes, 'FontSize', 18, 'LineWidth', 1.0, 'TickLength', [0.015 0.015], 'GridLineStyle', '-', ...
%      'XTick', [0 0.2 0.4 0.6 0.8 1],'YTick', [-0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3], 'YDir', 'Reverse');
% plot(x,Cf,'-b',x2,Cf2,'--r','LineWidth',1.5); 
% xlabel('$x/c$','FontSize', 24, 'Interpreter', 'latex');
% ylabel('$C_f$','FontSize', 24, 'Interpreter', 'latex');
% legend({'HDG', 'EDG'},'Location','SE','Interpreter','latex');
% box on;
% grid on;
% axis on;
% axis([0 1 -0.4 0.3]);
% print -depsc nstrefftzcfp2.eps
% 
