porder = 4;
ngrid = 21;
nstage = 2;
torder = 2;
hybrid = 'hdg';
clear app;

app.adjoint = 0;
app.denseblock = 0;
app.hybrid = hybrid;
app.localsolve=1;


tau = 1;
app.source = 'source2d';
app.flux = 'flux2d';
app.fbou = 'fbou';
app.fhat = 'fhat';
app.arg = {tau};
app.localsolve=1;

app.bcm = [3; 3; 3; 1;];
app.bcs = [0; 0; 0; 0];

app.bcd = [];
app.bcv = [];

app.hybrid = hybrid;
app.tdep = true;
app.wave = false;
app.alag = false;
app.flg_q = 1;
app.flg_p = 0;
app.flg_g = 0;

app.fc_q = 1;          
app.fc_u = 1;
app.fc_p = 0;


app.nd = 2;
app.nch  = 1;                       % Number of componets of UH
app.nc   = app.nch*(app.nd+1);    % Number of componeents of UDG
app.ncu = app.nch;

% dt = 1e-2;
ntime  = 100;
% dt = linspace(dt, dt*ntime, ntime);
dt = 0.01*ones(ntime,1);

app.time = [];
app.dtfc = [];
app.alpha = [];

% Initializing data structures




mesh = mkmesh_rect(41,81,porder,0,[0 1 0 2],0,1);
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);

% % Number density initialized to the same gaussian for both electrons and positives, and 0 for negatives.
% % These have to be initialized in the same order that UDG is in
% %                 ne_0,         q_ne_r0,        q_ne_z0
initu_func_set = {@initu_func; @initq_func_r; @initq_func_z;};
UDG = initu(mesh,initu_func_set,app.arg); 
UDG_history = zeros([size(UDG),ntime+1]);
UDG_history(:,:,:,1) = UDG;     % Adding in IC to the first snapshot
itime_restart = 0;
UH=inituhat(master,mesh.elcon,UDG,app.ncu);

% Plot IC
% u = UDG_history(:,1,:,1);
% scaplot(mesh,u,[],0,1); axis equal; axis tight; colormap jet; title('u0');
% return;





% porder = 4;
% ngrid  = 21;
% nstage = 2;
% torder = 2;
% hybrid = 'hdg';
% clear app;

% app.adjoint = 0;
% app.denseblock = 0;
% app.hybrid = hybrid;
% app.localsolve=1;
% ntime  = 50;
% dt = 0.01*ones(ntime,1);

% kappa = 1e-16;
% c = [0,1]; %1,8;2,5];
% tau = 1;

% app.source = 'source2d';
% app.flux = 'flux2d';
% app.fbou = 'fbou';
% app.fhat = 'fhat';

% app.arg = {kappa,c,tau};
% app.bcm = [3;3;3;1];
% app.bcs = [0;0;0;0]; %[1,0,0;1,0,0];
% app.bcd = [1;1;1;1];
% app.bcv = [0;0;0;0]; 

% app.tdep = true;
% app.wave = false;
% app.alag = false;
% app.flg_q = 1;
% app.flg_p = 0;
% app.flg_g = 0;
% app.fc_q = 1;
% app.fc_u = 1;

% app.np = 2;
% app.nd = 2;
% app.nch  = 1;                       % Number of componets of UH
% app.nc   = app.nch*(app.nd+1);    % Number of componeents of UDG
% app.ncu  = 1;

% app.time = [];
% app.dtfc = [];
% app.alpha = [];

% mesh   = mkmesh_rect((ngrid-1)/2+1,ngrid,porder,0,[0 1 -1 1],0,0);
% master = mkmaster(mesh,2*porder);
% [master,mesh] = preprocess(master,mesh,hybrid);

% x = mesh.dgnodes(:,1,:);
% y = mesh.dgnodes(:,2,:);
% r2 = x.^2 + y.^2;
% UDG = initu(mesh,{0;0;0});
% UDG(:,1,:) = exp(-50*r2);
% UDG(:,2,:) = 50*x.*exp(- 50*r2);
% UDG(:,3,:) = 50*y.*exp(- 50*r2);
% figure(1); clf; scaplot(mesh,UDG(:,1,:),[],2,1); axis equal; axis tight;  axis on;

% UH = inituhat(master,mesh.elcon,UDG,app.ncu);

time = 0;
disp('Starting sim...')

for itime = 1:ntime
    fprintf('Timestep :  %d\n', itime+itime_restart);

    [UDG,UH] = hdg_solve_dirk(master,mesh,app,UDG,UH,[],time,dt(itime),nstage,torder);
    time = time + dt(itime);

    % UDG_history(:,:,:,itime+1+itime_restart) = UDG;
    figure(1); clf; scaplot(mesh,UDG(:,1,:),[],2,1); axis equal; axis tight; axis on;  
    pause(0.2);
    
end
