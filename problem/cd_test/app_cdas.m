porder = 4;
hybrid = 'hdg';
nstage = 2; 
torder = 2;
tau = 1;
clear app;
app.axisymmetry = 1;

app.source = 'source2d';
app.flux = 'flux2d';
app.fbou = 'fbou_electrondensity';
app.fhat = 'fhat_electrondensity';
app.localsolve=1;

app.bcm = [2; 3; 2; 1;];
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

ntime  = 68;
dt = 0.01*ones(ntime,1);

app.time = [];
app.dtfc = [];
app.alpha = [];

% Initializing data structures
app.arg = {tau};

mesh = mkmesh_rect(41,81,porder,0,[0 1 0 2],0,1);
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);

% Number density initialized to the same gaussian for both electrons and positives, and 0 for negatives.
% These have to be initialized in the same order that UDG is in
%                 ne_0,         q_ne_r0,        q_ne_z0
initu_func_set = {@initu_func; @initq_func_r; @initq_func_z;};
UDG = initu(mesh,initu_func_set,app.arg); 

% Plot IC
% u = UDG_history(:,1,:,1);
% scaplot(mesh,u,[],0,1); axis equal; axis tight; colormap jet; title('u0');
% return;

UH=inituhat(master,mesh.elcon,UDG,app.ncu);

time = 0;
disp('Starting sim...')

for itime = 1:ntime
    fprintf('Timestep :  %d\n', itime+itime_restart);

    [UDG,UH] = hdg_solve_dirk(master,mesh,app,UDG,UH,[],time,dt(itime),nstage,torder);
    time = time + dt(itime);

    figure(1); clf; scaplot(mesh,UDG(:,1,:),[],2,1); axis equal; axis tight; axis on;  
    % pause(0.2);
    
end
