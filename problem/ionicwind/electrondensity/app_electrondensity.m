porder = 2;
hybrid = 'hdg';
nstage = 1;
torder = 1;
tau = 1;

app.source = 'source_electrondensity';
app.flux = 'flux_electrondensity';
app.fbou = 'fbou_electrondensity';
app.fhat = 'fhat_electrondensity';
app.localsolve=1;
% app.bcm = [4;4;4;2;2;2;1];
app.bcm = [3;3;3;1;4;4;3];
app.bcs = [0;0;0;0;0;0;0];
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
app.nch  = 1;                       % Number of componets of UH
app.nc   = app.nch*(app.nd+1);    % Number of componeents of UDG
app.ncu = 1;

ntime  = 40;
dt = linspace(.05, 2 ,ntime);

% Physical parameters
Kep = 2e-13;             % mu[1] Recombination coeff - pos and neg ions [m^3/s]
Knp = 2e-13;             % mu[2] Recombination coeff - pos ions and electrons [m^3/s]
mu_p = 2.43e-4;          % mu[3] Pos ion mobility [m^2/(Vs)]
mu_n = 2.7e-4;           % mu[4] Neg mobility [m^2/(Vs)]
De = 0.18;               % mu[5] Electron diffusion coefficient [m^2/s]
Dp = 0.028e-4;           % mu[6] Pos ion diffusion coefficient [m^2/s]
Dn = 0.043e-4;           % mu[7] Neg diffusion coefficient [m^2/s]
Nmax = 1e16;             % mu[8] Max number density for initial charge distribution [particles/m^3]
r0 = 0.0;                % mu[9] r-pos of emitter tip in reference frame [m]
z0 = 0.0;              % mu[10]z-pos of emitter tip in reference frame [m]
s0 = 2e-3;               % mu[11]Std deviation of initial charge distribution [m]
e = 1.6022e-19;          % mu[12]Charge on electron [C]
epsilon0 = 8.854e-12;     % mu[13]absolute permittivity of air [C^2/(N*m^2)]
Ua = -10e3;              % mu[14]Emitter potential relative to ground [V]
gamma = 0.001;           % mu[15]Secondary electron emission coefficient [1/m]
E_bd = 3e6;              % mu[16]Breakdown E field in air [V/m]
r_tip = 220e-6;          % mu[17] Tip radius of curvature [m]
D = 0.11736375131509072;     % Taken from swarm_params.py evaluated at E=3e6 V/m
n_ref = epsilon0*E_bd/(e*r_tip);  % Yes, this is redundant and could be recomputed from the above variables. But it saves having to recompute it each time in the functions.
% ^=7.5357e+17

P = 101325; % Pa
V = 1; % m^3
T = 273.15; % K
k_b = 1.380649e-23; % m2 kg s-2 K-1
N = P*V/(k_b*T);
mue_ref = 0.04266918357567234;   % m^2/(V*s)
D_star = D/(mue_ref*E_bd*r_tip);   % Nondimensionalized diffusion coefficient

% Set discretization parameters, physical parameters, and solver parameters
     %      1  2    3   4    5      6     7     8     9     10      11   12   13       14     15
app.arg = {r0, z0, s0, Nmax, e, epsilon0, Ua, gamma, E_bd, r_tip, n_ref, N, mue_ref, D_star, tau};

app.time = [];
app.dtfc = [];
app.alpha = [];

mesh   = mkmesh_chen(porder);
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);

load '../poissonsolution2.mat';
UDG = UDG*Ua/(E_bd*r_tip);      % Sign flip and scaling potential function for problem nondimensionalization
mesh.dgnodes(:,3,:) = UDG(:,2,:);
mesh.dgnodes(:,4,:) = UDG(:,3,:);

UDG = initu(mesh,{@initu_func;@initq_func_r;@initq_func_z},app.arg);
UDG_history = zeros([size(UDG),ntime+1]);
UDG_history(:,:,:,1) = UDG;
UH=inituhat(master,mesh.elcon,UDG,1);

time = 0;
figure(1);
% set(gcf,'color','black');
for itime = 1:ntime
    fprintf('Timestep :  %d\n', itime);
    
    [UDG,UH] = hdg_solve_dirk(master,mesh,app,UDG,UH,[],time,dt(itime),nstage,torder);
    time = time + dt(itime);
    UDG_history(:,:,:,itime+1) = UDG;
    
    % figure(1); scaplot(mesh,UDG(:,1,:),[],0,0); axis equal; axis tight; colormap jet; title('ne');
end

% [UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,0*UDG);      % CHANGE ME - IMPLEMENT TIMESTEPPING SCHEME

% clf;
% figure(1); scaplot(mesh,UDG(:,1,:),[],0,0); axis equal; axis tight; colormap jet; title('ne');
% figure(2); scaplot(mesh,UDG(:,2,:),[],0,0); axis equal; axis tight; colormap jet; title('dne/dr');
% figure(3); scaplot(mesh,UDG(:,3,:),[],0,0); axis equal; axis tight; colormap jet; title('dne/dz');

% figure(4); scaplot(mesh,mue*mesh.dgnodes(:,3,:),[],0,0); axis equal; axis tight; colormap jet; title('Er');
% figure(5); scaplot(mesh,mue*mesh.dgnodes(:,4,:),[],0,0); axis equal; axis tight; colormap jet; title('Ez');
