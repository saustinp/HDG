setapplicationpath('FM/condiff');

porder = 5;
ngrid  = 5;
hybrid = 'hdg';

kappa = 1;
c = [20,20]; 
tau = 20;

app.adjoint = 0;
app.denseblock = 0;
app.hybrid = 'hdg';
app.localsolve=1;
app.arg = {kappa,c,tau};
app.bcm = [1;1;1;1];
app.bcs = [0;0;0;0]; 
app.bcd = [1;1;1;1];
app.bcv = [0;0;0;0]; 

app.tdep = false;
app.wave = false;
app.alag = false;
app.flg_q = 1;
app.flg_p = 0;
app.flg_g = 0;

app.fc_q = 1;
app.fc_u = 0;
app.fc_p = 0;

app.np = 2;
app.nd = 2;
app.nch  = 1;                       % Number of componets of UH
app.nc   = app.nch*(app.nd+1);    % Number of componeents of UDG
app.ncu = 1;

app.time = [];
app.dtfc = [];
app.alpha = [];

mesh   = mkmesh_square(ngrid,ngrid,porder);
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);

UDG = initu(mesh,{0;0;0});
UH = inituhat(master,mesh.elcon,UDG,app.ncu);

% HDG solver
[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);
figure(1); clf; scaplot(mesh,UDG(:,1,:),[],2); axis equal; axis tight;

% % HDG postprocessing 
% mesh1 = mkmesh_square(ngrid,ngrid,porder+1);
% master1 = mkmaster(mesh1,2*(porder+1));
% [master1,mesh1] = preprocess(master1,mesh1,hybrid);
% UDGstar = postprocessnd(master,mesh,master1,mesh1,UDG);
% 
% %figure(2); clf; scaplot(mesh1,UDGstar(:,1,:),[],2); axis equal; axis tight;
% 
% [VDG,VH] = hdg_solve_adjoint(master,mesh,UDG,UH,[],app);
% figure(3);clf; scaplot(mesh,VDG(:,1,:),[],2); axis equal; axis tight;
% 

