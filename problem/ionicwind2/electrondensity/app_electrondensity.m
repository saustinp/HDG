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

dt = 1e-3;
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
E_bd = app.arg{15};
r_tip = app.arg{16};
UDG = UDG*Ua/(E_bd*r_tip);      % Sign flip and scaling potential function for problem nondimensionalization
electrostatic_sol = UDG;
mesh.dgnodes(:,3,:) = UDG(:,2,:);
mesh.dgnodes(:,4,:) = UDG(:,3,:);

% scaplot(mesh,electrostatic_sol(:,1,:),[-15.15 0],0,0); title('Phi');
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

UH=inituhat(master,mesh.elcon,UDG,app.ncu);

time = 0;
disp('Starting sim...')

for itime = 1:ntime
    fprintf('Timestep :  %d\n', itime);
    
    [UDG,UH] = hdg_solve_dirk(master,mesh,app,UDG,UH,[],time,dt(itime),nstage,torder);
    time = time + dt(itime);

    % UDG(:,1:3,:) = max(UDG(:,1:3,:), 0);     % Clip to 0 to prevent wiggles

    UDG_history(:,:,:,itime+1) = UDG;
    UDG_history(:,4,:,itime+1) = UDG_history(:,4,:,itime+1) + electrostatic_sol(:,1,:);
    UDG_history(:,8,:,itime+1) = UDG_history(:,8,:,itime+1) + electrostatic_sol(:,2,:);
    UDG_history(:,12,:,itime+1) = UDG_history(:,12,:,itime+1) + electrostatic_sol(:,3,:);
    

    % % Plotting the E field
    % Er0 = mesh.dgnodes(:,3,:);
    % Ez0 = mesh.dgnodes(:,4,:);
    % Er_prime = UDG_history(:,8,:,itime+1);
    % Ez_prime = UDG_history(:,12,:,itime+1);
    % Er = Er_prime + Er0;
    % Ez = Ez_prime + Ez0;
    % normE = sqrt(Er.^2 + Ez.^2);
    % N = 2.4614924955148245e25;

    % scaplot(mesh,normE/N/1e-21,[],0,0); axis equal; axis tight; colormap jet; title('');
    % return;
end

% animate_sol;

% [UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,0*UDG);      % CHANGE ME - IMPLEMENT TIMESTEPPING SCHEME

% clf;
% figure(1); scaplot(mesh,UDG(:,1,:),[],0,0); axis equal; axis tight; colormap jet; title('ne');
% figure(2); scaplot(mesh,UDG(:,2,:),[],0,0); axis equal; axis tight; colormap jet; title('dne/dr');
% figure(3); scaplot(mesh,UDG(:,3,:),[],0,0); axis equal; axis tight; colormap jet; title('dne/dz');

% figure(4); scaplot(mesh,mue*mesh.dgnodes(:,3,:),[],0,0); axis equal; axis tight; colormap jet; title('Er');
% figure(5); scaplot(mesh,mue*mesh.dgnodes(:,4,:),[],0,0); axis equal; axis tight; colormap jet; title('Ez');
