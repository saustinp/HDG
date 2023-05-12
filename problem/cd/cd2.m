setapplicationpath('FM/condiff');

porder = 4;
ngrid  = 9;
nstage = 1;
torder = 1;
hybrid = 'hdg';

app.adjoint = 0;
app.denseblock = 0;
app.hybrid = hybrid;
app.localsolve=1;
ntime  = 5;
dt = linspace(0.01,10,ntime);

kappa = 0.5;
c = [10,10]; %1,8;2,5];
tau = 10;

app.arg = {kappa,c,tau};
app.bcm = [1;1;1;1];
app.bcs = [0;0;0;0]; %[1,0,0;1,0,0];
app.bcd = [1;1;1;1];
app.bcv = [0;0;0;0]; 

app.tdep = true;
app.wave = false;
app.alag = false;
app.flg_q = 1;
app.flg_p = 0;
app.flg_g = 0;

app.np = 2;
app.nd = 2;
app.nch  = 1;                       % Number of componets of UH
app.nc   = app.nch*(app.nd+1);    % Number of componeents of UDG
app.ncu  = 1;

app.time = [];
app.dtfc = [];
app.alpha = [];

mesh   = mkmesh_square(ngrid,ngrid,porder);
master = mkmaster(mesh,2*porder);

[master,mesh,app] = preprocess(master,mesh,app);

UDG = initu(mesh,{0;0;0});
UH = inituhat(master,mesh.elcon,UDG,app.ncu);

time = 0;
figure(1);
set(gcf,'color','black');
for itime = 1:ntime
    fprintf('Timestep :  %d\n', itime);
    
    [UDG,UH] = hdg_solve_dirk(master,mesh,app,UDG,UH,[],time,dt(itime),nstage,torder);    
    time = time + dt(itime); 
    
    figure(1); scaplot(mesh,UDG(:,1,:),[],2,1); axis equal; axis tight;    
end

% HDG solver
app.fc_q = 1;
app.fc_p = 0;
app.tdep = false;
[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);

% HDG postprocessing 
mesh1 = mkmesh_square(ngrid,ngrid,porder+1);
master1 = mkmaster(mesh1,2*(porder+1));
[master1,mesh1] = preprocess(master1,mesh1,hybrid);
UDGstar = postprocessnd(master,mesh,master1,mesh1,UDG);

figure(1); scaplot(mesh,UDG(:,1,:),[],2,1); axis equal; axis tight;
figure(2); scaplot(mesh1,UDGstar(:,1,:),[],2,1); axis equal; axis tight;
