setapplicationpath('FM/euler')

hybrid = 'hdg';

m      = 20;
n      = 40;
gridNum = 1;
porder = 3;
torder = 1;
nstage = 1;

gam = 1.4;
epslm = 0.00;
Minf = 0.4;                  % Infinity conditions
pinf = 1/(gam*Minf^2);
alpha = 2*(pi/180);
Re = 1000;  % Irrelevant value
Pr = 0.73; % Irrelevant value
tau = 2;
ui = [1, cos(alpha), sin(alpha), 0.5+pinf/(gam-1)];

app.source = 'source';
app.flux = 'flux';
app.fbou = 'fbou';
app.fhat = 'fhat';
app.iterative=0;
app.hybrid = hybrid;
app.localsolve=0;
app.arg = {gam,epslm,Re,Pr,Minf,tau};

% app.bcm  = [7,1];  % 2: Slip wall, 1: Far-field
app.bcm  = [2,1];  % 2: Slip wall, 1: Far-field
app.bcs  = [ui; ui];

app.bcd  = [1,3];  % 2: Slip wall, 1: Far-field
app.bcv  = [ui; ui];

app.wave = false;
app.tdep = true;
app.alag = false;
app.flg_q = 0;
app.flg_p = 0;
app.flg_g = 0;

app.np = 2;
app.nd = 2;
app.nch = 2+app.nd;                % Number of componets of UH
app.ncq = 0;
app.ncu = app.nch;                        % Number od components with time derivative
app.nc  = app.ncu;                 % Number of componeents of UDG
app.ncp = 0;

app.dtfc = 0;
app.alpha = 0;

app.adjoint = 0;
app.linearproblem = 0;
app.appname = 'euler';
app.linearSolver = 1;
app.jacobianStep = 0;
app.orderingStep = 0;
app.dt = [1e-4,1e-3,1e-2,0.1,1,10,100];
app.torder = torder;
app.nstage = nstage;
app.time = 0;
app.fc_q = 0;
app.fc_u = 0;
app.fc_p = 0;
app.ns   = 1;

% mesh = mkmesh_naca0012(porder,1,gridNum);
mesh = mkmesh_trefftz(m,n,porder,[0.05,0.05,1.98]);
master = mkmaster(mesh,2*porder);
[master,mesh,app] = preprocess(master,mesh,app);

UDG0 = initu(mesh,{ui(1),ui(2),ui(3),ui(4)});
UH0 = inituhat(master,mesh.elcon,UDG0,app.ncu);

fprintf('\nSteady State Solve\n');
[UDG,UH] = hdg_solve(master,mesh,app,UDG0,UH0,[]);

figure(1); clf;
scaplot(mesh,eulereval(UDG(:,1:4,:),'M',1.4),[],3,0); 
hold on;
axis off; axis equal; axis tight;    
colormap('jet');
colorbar('FontSize',14);
axis([-0.2 1.2 -0.6 0.6]);


return;

% %DRIVER FOR THE EULER EQUATIONS IN A SUBSONIC FLOW OVER TREFFTZ FOIL
% 
% setapplicationpath('FM/euler');
% 
% m      = 20;
% n      = 40;
% porder = 4;
% nstage = 1;
% torder = 1;
% hybrid = 'hdg';
% 
% gam = 1.4;
% epslm = 0.0;
% tau = 2;
% rmin = 0.1;
% Minf = 0.2;                  % Infinity conditions
% pinf = 1/(gam*Minf^2);
% alpha = 0*pi/180;
% ui = [ 1, cos(alpha), sin(alpha), 0.5+pinf/(gam-1)];
% 
% % ntime = 10;
% % dt = linspace(0.01,2,ntime);
% dt = [0.01 0.1 1 10];
% ntime = length(dt);
% 
% app.source = 'source';
% app.flux = 'flux';
% app.fbou = 'fbou';
% app.fhat = 'fhat';
% app.iterative=0;
% app.hybrid = hybrid;
% app.localsolve=0;
% app.arg = {gam,epslm,tau};
% app.bcm  = [2,1];  % 2: Wall, 1: Far-field, Bottom-Top: Wall, Left-Right: Far-field
% app.bcs  = [ui; ui];
% 
% app.bcd  = [1,3];  % 2: Wall, 1: Far-field
% app.bcv  = [ui; ui];
% 
% app.wave = false;
% app.tdep = true;
% app.alag = false;
% app.flg_q = 0;
% app.flg_p = 0;
% app.flg_g = 0;
% 
% app.np = 2;
% app.nd = 2;
% app.nch  = 2+app.nd;                % Number of componets of UH
% app.nc   = app.nch;                 % Number of componeents of UDG
% app.ncu  = 4;                        % Number od components with time derivative
% 
% app.time = [];
% app.dtfc = [];
% app.alpha = [];
% 
% app.adjoint = 0;
% app.linearproblem = 0;
% app.appname = 'euler';
% app.linearSolver = 1;
% app.jacobianStep = 0;
% app.orderingStep = 0;
% app.dt = dt;
% app.torder = torder;
% app.nstage = nstage;
% app.time = 0;
% app.fc_q = 0;
% app.fc_u = 0;
% app.fc_p = 0;
% app.ns   = 1;
% 
% mesh = mkmesh_trefftz(m,n,porder,[0.05,0.05,1.98]);
% master = mkmaster(mesh,2*porder);
% [master,mesh,app] = preprocess(master,mesh,app);
% 
% UDG = initu(mesh,{ui(1),ui(2),ui(3),ui(4)});
% UH = inituhat(master,mesh.elcon,UDG,app.ncu);
% 
% tic
% time = 0;
% figure(1);
% set(gcf,'color','black');
% for itime = 1:ntime
%     fprintf('Timestep :  %d\n', itime);
% 
%     %[Un,Hn,Pn] = hdg_solve_dirk(master,mesh,app,UDG,UH,PDG,time,dt,nstage,torder)
%     [UDG,UH] = hdg_solve_dirk(master,mesh,app,UDG,UH,[],time,dt(itime),nstage,torder);
%     time = time + dt(itime); 
%     
%     scaplot(mesh,eulereval(UDG(:,1:app.nch,:),'M',gam),[],2,1); axis off; %axis equal; axis tight;
%     tm=sprintf('t* = %.3f',time);
%     title(tm,'Color','white','FontSize',16,'FontName','Courier');
% %     set(gcf,'InvertHardCopy','off');   % conserve the black background when printing image
% %     pname=sprintf('plot/euler_%05d.png',itime);  % file name
% %     screen2png(pname)
% end
% 
% fprintf('\nSteady State Solve\n');
% app.fc_q = 1;
% app.fc_u = 0;
% app.tdep = false;
% app.adjoint = 0;
% [UDG,UH]=hdg_solve(master,mesh,app,UDG,UH,[]);
% 
% figure(1); clf; scaplot(mesh,eulereval(UDG(:,1:app.nch,:),'M',gam),[],2); axis off;
% 
% app.adjoint = 1;
% [VDG,VH] = hdg_solve(master,mesh,app,UDG,UH,[]);
% figure(2); clf; scaplot(mesh,VDG(:,1,:),[],2); axis off;



% save trefftzedg.mat master mesh app UDG UH
% 
% load trefftzedg.mat
% app.arg={1.4,0,0,0,Minf};
% [Cp1,Cf1,x1]=getsurfacedata(master,mesh,app,UDG,UH,2);
% 
% load trefftzhdg.mat
% app.arg={1.4,0,0,0,Minf};
% [Cp2,Cf2,x2]=getsurfacedata(master,mesh,app,UDG,UH,2);
% 
% figure(1); clf;
% set(axes, 'FontSize', 18, 'LineWidth', 1.0, 'TickLength', [0.015 0.015], 'GridLineStyle', '-', ...
%      'XTick', [0 0.2 0.4 0.6 0.8 1],'YTick', [-1.0:0.2:0.4], 'YDir', 'Reverse');
% plot(x1,Cp1,'-b',x2,Cp2,'--r','LineWidth',1.5); 
% xlabel('$x/c$','FontSize', 24, 'Interpreter', 'latex');
% ylabel('$-C_p$','FontSize', 24, 'Interpreter', 'latex');
% legend({'EDG', 'HDG'},'Location','SE','Interpreter','latex');
% box on;
% grid on;
% axis on;
% axis([0 1 -1 0.5]);
% print -depsc eulertrefftcp.eps
% 
% 
% 
% 
% 
% 
% 
