porder = 2;
hybrid = 'hdg';
kappa = 0.1;
c = [1, -1];
tau = 1;

app.source = 'source_cd';
app.flux = 'flux_cd';
app.fbou = 'fbou_cd';
app.fhat = 'fhat_cd';
app.localsolve=1;
app.arg = {kappa,c,tau};
app.bcm = [4;2;1;2;2;2;1];
app.bcs = [0;0;1;0;0;0;1];
app.bcd = [];
app.bcv = [];

app.hybrid = hybrid;
app.tdep = false;
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

app.time = [];
app.dtfc = [];
app.alpha = [];

mesh   = mkmesh_chen(porder);
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);

UDG = initu(mesh,{0;0;0});
UH=inituhat(master,mesh.elcon,UDG,1);

[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,0*UDG);

figure(1); scaplot(mesh,UDG(:,1,:),[],0,1); axis equal; axis tight; colormap jet;
figure(2); scaplot(mesh,UDG(:,2,:),[],0,1); axis equal; axis tight; colormap jet;
figure(3); scaplot(mesh,UDG(:,3,:),[],0,1); axis equal; axis tight; colormap jet;







