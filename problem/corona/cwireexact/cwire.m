setapplicationpath('FM/corona');

porder = 4;
ngrid  = 15;

K = 1;
D = 0;
e0 = 1;
Ec = sqrt(2);
a = 1;
b = 10;
J0 = 2*pi;
tau = 1;
param = {K,D,e0,Ec,a,b,J0,tau};

hybrid = 'hdg';
app.source = 'source';
app.flux = 'flux';
app.fbou = 'fbou';
app.fhat = 'fhat';
app.adjoint = 0;
app.denseblock = 0;
app.hybrid = hybrid;
app.localsolve=1;
app.arg = param;
app.bcm = [1;3];
app.bcs = [0;0]; 
app.bcd = [1;1];
app.bcv = [0;0]; 

app.denseblock = 0;
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
app.nch  = 2;                       % Number of componets of UH
app.nc   = app.nch*(app.nd+1);    % Number of componeents of UDG
app.ncu = 2;

app.appname = 'corona';
app.time = [];
app.dtfc = [];
app.alpha = [];

mesh = mkmesh_circincirc(porder,ngrid,ngrid,a,b);
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,app);
%return;

UDG = initu(mesh,{1;1;0;0;0;0});
UH = inituhat(master,mesh.elcon,UDG,app.ncu);

% HDG solver
[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);

for i = 1:6
    figure(i); clf; scaplot(mesh,UDG(:,i,:),[],2); 
    axis equal; axis tight; axis off; colormap jet;
    colorbar('FontSize',15);    
%     fn = ['numsol' num2str(i)];
%     print('-dpng',fn);
end

UEX = exactsol(mesh.dgnodes, cell2mat(param));
for i = 1:6
    figure(6+i); clf; scaplot(mesh,abs(UEX(:,i,:)-UDG(:,i,:)),[],2); 
    axis equal; axis tight; axis off; colormap jet;
    colorbar('FontSize',15);
%     fn = ['errsol' num2str(i)];
%     print('-dpng',fn);
end








