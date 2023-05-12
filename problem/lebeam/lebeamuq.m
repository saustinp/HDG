%DRIVER FOR CONVECTION DIFFUSION
%
setapplicationpath('SM/LE_uq');

porder = 3;
ngrid  = 5;

distortion = 0.0;

mu = 1;
lambda = [10]; %1,8;2,5];
tau = 1;
alpha = 0.01;

hybrid = 'hdg';
app.localsolve=1;
app.arg = {mu,lambda,tau};
app.bcm = [3;3;3;1];
app.bcs = [0,0,0;0,-0.5,0;0,0,0;0,0,0]; %[1,0,0;1,0,0];
app.bcd = [1;1;1;1];
app.bcv = [0;0;0;0]; 

app.denseblock = 0;
app.tdep = false;
app.wave = false;
app.alag = false;
app.flg_q = 1;
app.flg_p = 1;
app.flg_g = 0;

app.fc_q = 1;
app.fc_u = 0;
app.fc_p = 1;

app.np = 2;
app.nd = 2;
app.nch  = 2;                       % Number of componets of UH
app.nc   = app.nch*(app.nd+1)+1;    % Number of componeents of UDG
app.ncu = 2;

app.time = [];
app.dtfc = [];
app.alpha = alpha;

mesh = mkmesh_square(ngrid,ngrid,porder,0,2,1);
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);

UDG = initu(mesh,{0;0;0;0;0;0});
UH = inituhat(master,mesh.elcon,UDG,app.ncu);

% HDG solver
tic
[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);
toc

mesh1=mesh;
mesh1.dgnodes(:,1,:)=mesh.dgnodes(:,1,:)+UDG(:,1,:);
mesh1.dgnodes(:,2,:)=mesh.dgnodes(:,2,:)+UDG(:,2,:);

figure(1); scaplot(mesh,real(UDG(:,1,:)),[],2,1); axis equal; axis tight;
figure(2); scaplot(mesh,real(UDG(:,3,:)),[],2,1); axis equal; axis tight;
figure(3); meshplot(mesh1,1); axis equal; axis tight;



