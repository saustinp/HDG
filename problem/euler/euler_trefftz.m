%DRIVER FOR THE EULER EQUATIONS IN A SUBSONIC FLOW OVER TREFFTZ FOIL

setapplicationpath('FM/euler');

m      = 15;
n      = 30;
porder = 3;
nstage = 1;
torder = 1;
hybrid = 'hdg';

gam = 1.4;
epslm = 0.0;
tau = 2;
rmin = 0.1;
Minf = 0.2;                  % Infinity conditions
pinf = 1/(gam*Minf^2);
alpha = 0*pi/180;
ui = [ 1, cos(alpha), sin(alpha), 0.5+pinf/(gam-1)];

% ntime = 10;
% dt = linspace(0.01,2,ntime);
dt = [0.01 0.1 1 10];
ntime = length(dt);

clear app;
app.iterative=0;
app.hybrid = hybrid;
app.localsolve=0;
app.arg = {gam,epslm,tau,rmin};
app.bcm  = [2,1];  % 2: Wall, 1: Far-field, Bottom-Top: Wall, Left-Right: Far-field
app.bcs  = [ui; ui];

app.bcd  = [1,3];  % 2: Wall, 1: Far-field
app.bcv  = [ui; ui];

app.wave = false;
app.tdep = true;
app.alag = false;
app.flg_q = 0;
app.flg_p = 0;
app.flg_g = 0;

app.np = 2;
app.nd = 2;
app.nch  = 2+app.nd;                % Number of componets of UH
app.nc   = app.nch;                 % Number of componeents of UDG
app.ncu  = 4;                        % Number od components with time derivative

app.time = [];
app.dtfc = [];
app.alpha = [];

mesh = mkmesh_trefftz(m,n,porder,[0.05,0.05,1.98]);
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);

UDG = initu(mesh,{ui(1),ui(2),ui(3),ui(4)});
UH = inituhat(master,mesh.elcon,UDG,app.ncu);


fprintf('\nSteady State Solve\n');
app.fc_q = 1;
app.fc_u = 0;
app.tdep = false;
app.adjoint = 0;
[UDG,UH]=hdg_solve(master,mesh,app,UDG,UH,[]);

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
% 
% figure(1); clf; scaplot(mesh,eulereval(UDG(:,1:app.nch,:),'M',gam),[],2); axis off;
% 
% app.adjoint = 1;
% [VDG,VH] = hdg_solve(master,mesh,app,UDG,UH,[]);
% figure(2); clf; scaplot(mesh,VDG(:,1,:),[],2); axis off;


