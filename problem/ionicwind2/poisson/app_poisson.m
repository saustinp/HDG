porder = 2;
hybrid = 'hdg';
kappa = 1;
tau = 1;

app.source = 'source_poisson';
app.flux = 'flux_poisson';
app.fbou = 'fbou_poisson';
app.fhat = 'fhat_poisson';
app.localsolve=1;
app.arg = {kappa,tau};
app.bcm = [3;3;3;1;1;1;1];
app.bcs = [0;0;0;0;0;0;1];
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
app.nc   = app.nch*(app.nd+1);    % Number of components of UDG: for each equation, one solution variable and ND number of gradient variables
app.ncu = 1;

app.time = [];
app.dtfc = [];
app.alpha = [];

mesh   = mkmesh_chen(porder, "chen_19k.msh");
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);

UDG = initu(mesh,{0;0;0});
UH=inituhat(master,mesh.elcon,UDG,1);

[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,0*UDG);

save '../poissonsolution_19k.mat' UDG

clf;
% figure(1); scaplot(mesh,UDG(:,1,:),[],0,0); axis equal; axis tight; colormap turbo; title('Phi');
% figure(2); scaplot(mesh,UDG(:,2,:),[],0,0); axis equal; axis tight; colormap turbo; title('Er');
% figure(3); scaplot(mesh,UDG(:,3,:),[],0,0); axis equal; axis tight; colormap turbo; title('Ez');

normE = sqrt(UDG(:,2,:).^2 + UDG(:,3,:).^2);
figure(); scaplot(mesh,normE,[],0,0); axis equal; axis tight; colormap turbo; title('||E||, Vapp=1V');
N=2.4614924955148245e25;
Vapp=10e3;
EN_Td = Vapp*normE/N/1e-21;
figure(); scaplot(mesh,EN_Td,[],0,0); axis equal; axis tight; colormap turbo; title('E/N [Td], Vapp=10kV');

EN_Td_thresh = (EN_Td>130)+1;
figure(); scaplot(mesh,EN_Td_thresh,[],0,0); axis equal; axis tight; colormap turbo; title('E/N [Td] > 130, mask');