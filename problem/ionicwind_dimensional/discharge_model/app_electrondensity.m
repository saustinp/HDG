porder = 2;
hybrid = 'hdg';
nstage = 2;
torder = 2;
tau = 1;
app.axisymmetry = 1;

app.source = 'source2d';
app.flux = 'flux2d';
app.fbou = 'fbou_electrondensity';
app.fhat = 'fhat_electrondensity';
app.localsolve=1;

% Boundaries
% 1 Axisymmetry
% 2 Right farfield
% 3 Top farfield
% 4 Outflow
% 5 Ground plane
% 6 Cylinder
% 7 Needle

% BC types
% 1 Dirichlet scalar
% 3 Prescribed flux - scalar
% 4 Diffusion flux only
% 5 Symmetry
% 6 Prescribed flux - function

% Leaving out the electrostatic equation
app.bcm = [1; 2; 2; 3; 4; 4; 5];
app.bcs = [0; 0; 0; 0; 0; 0; 0];

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
app.fc_u = 0;
app.fc_p = 0;

app.nd = 2;
app.nch  = 4;                       % Number of componets of UH
app.nc   = app.nch*(app.nd+1);    % Number of componeents of UDG
app.ncu = app.nch;

dt = 1e-8;
ntime  = 20;
dt = linspace(dt, dt*ntime, ntime);

app.time = [];
app.dtfc = [];
app.alpha = [];

% Initializing data structures
app.arg = init_phys_param();     % Physics param loaded in a separate script
app.arg{end+1} = tau;

mesh = mkmesh_chen(porder, "chen_11k.msh");
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);

load '../poissonsolution_11k.mat';
Ua = app.arg{13};
UDG = UDG*Ua;
electrostatic_sol = UDG;
mesh.dgnodes(:,3,:) = UDG(:,2,:);
mesh.dgnodes(:,4,:) = UDG(:,3,:);

% return;
% scaplot(mesh,electrostatic_sol(:,1,:),[],0,0); title('Phi');
% return;

% Number density initialized to the same gaussian for both electrons and positives, and 0 for negatives.
% These have to be initialized in the same order that UDG is in
%                 ne_0, nn_0, np_0, phi_0,      q_ne_r0, q_nn_r0, q_np_r0, q_phi_r0,   q_ne_z0, q_nn_z0, q_np_z0, q_phi_z0
initu_func_set = {@initu_func;0;@initu_func;0;    @initq_func_r;0;@initq_func_r;0;    @initq_func_z;0;@initq_func_z;0};
UDG = initu(mesh,initu_func_set,app.arg); 

% Creating history vector to track the snapshots. Set the first snapshot equal to the initial condition as a test
UDG_history = zeros([size(UDG),ntime+1]);
UDG_history(:,:,:,1) = UDG;     
UDG_history(:,:,:,1) = UDG;
UDG_history(:,4,:,1) = UDG_history(:,4,:,1) + electrostatic_sol(:,1,:);
UDG_history(:,8,:,1) = UDG_history(:,8,:,1) + electrostatic_sol(:,2,:);
UDG_history(:,12,:,1) = UDG_history(:,12,:,1) + electrostatic_sol(:,3,:);

% return;

UH=inituhat(master,mesh.elcon,UDG,app.ncu);

time = 0;
disp('Starting sim...')

for itime = 1:ntime
    fprintf('Timestep :  %d\n', itime);
    
    [UDG,UH] = hdg_solve_dirk(master,mesh,app,UDG,UH,[],time,dt(itime),nstage,torder);
    time = time + dt(itime);

    UDG_history(:,:,:,itime+1) = UDG;
    UDG_history(:,4,:,itime+1) = UDG_history(:,4,:,itime+1) + electrostatic_sol(:,1,:);
    UDG_history(:,8,:,itime+1) = UDG_history(:,8,:,itime+1) + electrostatic_sol(:,2,:);
    UDG_history(:,12,:,itime+1) = UDG_history(:,12,:,itime+1) + electrostatic_sol(:,3,:);
end