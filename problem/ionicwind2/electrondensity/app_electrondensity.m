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

dt = 1e-4;
ntime  = 1000;
dt = linspace(dt, dt*ntime, ntime);

app.time = [];
app.dtfc = [];
app.alpha = [];

% Initializing data structures
app.arg = init_phys_param();     % Physics param loaded in a separate script
app.arg{end+1} = tau;

mesh = mkmesh_chen(porder, "chen_19k.msh");
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);

load '../poissonsolution_19k.mat';
Ua = app.arg{13};
E_bd = app.arg{15};
r_tip = app.arg{16};
UDG = UDG*Ua/(E_bd*r_tip);      % Sign flip and scaling potential function for problem nondimensionalization. Because the mesh is already nondimensionalized, this scaling factor works for both the potential and the E field
mesh.dgnodes(:,3,:) = UDG(:,2,:);
mesh.dgnodes(:,4,:) = UDG(:,3,:);

% Number density initialized to the same gaussian for both electrons and positives, and 0 for negatives.
% These have to be initialized in the same order that UDG is in
%                 ne_0, nn_0, np_0, phi_0,      q_ne_r0, q_nn_r0, q_np_r0, q_phi_r0,   q_ne_z0, q_nn_z0, q_np_z0, q_phi_z0
initu_func_set = {@initu_func;0;@initu_func;0;    @initq_func_r;0;@initq_func_r;0;    @initq_func_z;0;@initq_func_z;0};
UDG = initu(mesh,initu_func_set,app.arg); 

% Creating history vector to track the snapshots. Set the first snapshot equal to the initial condition as a test
UDG_history = zeros([size(UDG),ntime+1]);
UDG_history(:,:,:,1) = UDG;     % Adding in IC to the first snapshot

UH=inituhat(master,mesh.elcon,UDG,app.ncu);

time = 0;
disp('Starting sim...')

for itime = 1:ntime
    fprintf('Timestep :  %d\n', itime);

    % [stot, se, sn, sp, sphi] = source_eval(mesh.dgnodes,UDG,app.arg,0);
    % [f, cr_e, cz_e, cr_n, cz_n, cr_p, cz_p, De_s, Dn_s, Dp_s] = flux_eval(mesh.dgnodes,UDG,app.arg,0);
    % 
    % r1 = mesh.dgnodes(:,1,:);

    % Plot the electron convective velocities
    % figure();
    % scaplot(mesh, cr_e,[],0,0); axis equal; axis tight; colormap jet; title('cr_e');    
    % figure();
    % scaplot(mesh, cz_e,[],0,0); axis equal; axis tight; colormap jet; title('cz_e');
    % figure();
    % scaplot(mesh, cr_n,[],0,0); axis equal; axis tight; colormap jet; title('cr_n');    
    % figure();
    % scaplot(mesh, cz_n,[],0,0); axis equal; axis tight; colormap jet; title('cz_n');
    % figure();
    % scaplot(mesh, cr_p,[],0,0); axis equal; axis tight; colormap jet; title('cr_p');    
    % figure();
    % scaplot(mesh, cz_p,[],0,0); axis equal; axis tight; colormap jet; title('cz_p');

    % Plot the diffusion fields
    % figure();
    % scaplot(mesh, De_s,[],0,0); axis equal; axis tight; colormap jet; title('De_s');

    % Source term
    % figure();
    % scaplot(mesh, se,[],0,0); axis equal; axis tight; colormap jet; title('se');
    % figure();
    % scaplot(mesh, sn,[],0,0); axis equal; axis tight; colormap jet; title('sn');
    % figure();
    % scaplot(mesh, sp,[],0,0); axis equal; axis tight; colormap jet; title('sp');
    % figure();
    % scaplot(mesh, sphi,[],0,0); axis equal; axis tight; colormap jet; title('sphi');


    % Er = UDG_history(:,8,:,itime+1) + mesh.dgnodes(:,3,:);
    % Ez = UDG_history(:,12,:,itime+1) + mesh.dgnodes(:,4,:);
    % 
    % normE = sqrt(Er.^2+Ez.^2);
    % alpha = get_alpha(normE*E_bd, N);
    % eta = get_eta(normE*E_bd, N);
    % mue = get_mue(normE*E_bd, N);
    % diff = get_diffusion_e(normE*E_bd, N);

    % Plotting swarm params
    % figure();
    % scaplot(mesh, alpha,[],0,0); axis equal; axis tight; colormap jet; title('alpha');
    % figure();
    % scaplot(mesh, eta,[],0,0); axis equal; axis tight; colormap jet; title('eta');
    % figure();
    % scaplot(mesh, mue,[],0,0); axis equal; axis tight; colormap jet; title('mue');
    % figure();
    % scaplot(mesh, diff,[],0,0); axis equal; axis tight; colormap jet; title('diff');

    [UDG,UH] = hdg_solve_dirk(master,mesh,app,UDG,UH,[],time,dt(itime),nstage,torder);
    time = time + dt(itime);

    UDG_history(:,:,:,itime+1) = UDG;
    
end
