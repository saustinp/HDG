%setapplicationpath('FM/burgersav');

porder = 4;
nstage = 1;
torder = 1;
ngrid = 12;
hybrid = 'hdg';

clear app;
app.source = 'source';
app.flux = 'flux';
app.fbou = 'fbou';
app.fhat = 'fhat';
app.hybrid = hybrid;
app.localsolve=1;
ntime  = 5;
dt = linspace(0.01,10,ntime);

kappa = 0.0;
c = [0.5,0.0];
tau = 0.5;
av = 0.0;

app.arg = {kappa,c,av,1/(ngrid*porder),tau};
app.bcm = [4;1;2;1];
app.bcs = [0;0;0;0]; %[1,0,0;1,0,0];
app.bcd  = [1,1,1,1];  % 2: Slip wall, 1: Far-field
app.bcv  = [0; 0; 0; 0];

app.tdep = false;
app.wave = false;
app.alag = false;
app.flg_q = 1;
app.flg_p = 0;
app.flg_g = 0;

app.nd = 2;
app.nch  = 1;                       % Number of componets of UH
app.nc   = app.nch*(app.nd+1);    % Number of componeents of UDG
app.ncu  = 1;

app.time = [];
app.dtfc = [];
app.alpha = [];
app.itmax = 1;
app.fc_q = 1;
app.fc_p = 0;
app.fc_u = 0;
app.tdep = false;
app.nc = 3;

% mesh and master
mesh = mkmesh_shock2(porder);
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);
mesh = mkcgmesh(mesh);
mesh.ib = [];
mesh.in = 1:size(mesh.p2,1);
[ne, npe] = size(mesh.t2);

% initial solution
UDG = initu(mesh,{0.0;0;0});
UH=inituhat(master,mesh.elcon,UDG,1);

% HDG solver for constant viscosity field
mesh.dgnodes(:,3,:) = 0.01;
[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);
figure(1); clf;scaplot(mesh,UDG(:,1,:),[],2,0); axis equal; axis tight;

A = tensorproduct(master.plocvl,porder);
[cgnodes,cgelcon,rowent2elem,colent2elem,cgent2dgent] = mkcgent2dgent(mesh.dgnodes(:,1:2,:),1e-8);

beta = 0;
alpha = 100;
av = [1 1/2 1/4 1/8 1/16 1/32]*0.0025;
for i = 1:length(av)
    div = app.arg{4}*(UDG(:,2,:) + beta*UDG(:,3,:));
    s = ucg2udg(udg2ucg(div, cgent2dgent, rowent2elem), cgelcon);
    s = cgpoisson(mesh,master,s,[0.000001 1.0]);    
    s = s(mesh.t2');
    s = reshape(s,npe,1,ne);
    a = s.*(atan(alpha*s)/pi + 0.5) - atan(alpha)/pi + 0.5;
    mesh.dgnodes(:,3,:) = av(i)*a;    
            
    UDGprev = UDG;
    UHprev = UH;    
    [UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);
        
    [~,~,e1,e2] = dgprojection2(master,mesh,UDG(:,1,:),porder-1);
    [max(e1) max(e2)]
    
    u=A\squeeze(UDG(:,1,:));
    u1=u(1:(porder)*(porder),:);
    e1=sum(abs(u1),1);
    e=sum(abs(u),1);
    idx = e>1e-6;
    
%     [i av(i) min(e1(idx)./e(idx))]       
%     if min(e1(idx)./e(idx)) < 0.5
%         break;
%     end
end

figure(1); clf;scaplot(mesh, mesh.dgnodes(:,3,:) , [],2,0); axis equal; axis tight;
figure(2); clf;scaplot(mesh,UDG(:,1,:),[],2,0); axis equal; axis tight;       
figure(3); clf;scaplot(mesh,UDGprev(:,1,:),[],2,0); axis equal; axis tight;       

x = linspace(-1, 1, 1000);
y = 0.1:0.1:0.9;
ux = zeros(length(x), length(y));
for i = 1:length(y)
    xy = [x(:) y(i)*ones(length(x),1)];
    ux(:,i) = fieldatx(mesh,UDGprev(:,1,:),xy,20);
end
figure(4); clf; plot(x, ux, '-');

alpha = 2; beta = 4; gamma = 2;
ue = zeros(length(x), length(y));
for i = 1:length(y)        
    [up, um] = exactsolution(x, y(i), alpha, beta, gamma);
    ue(:,i) = um;
end

figure(4); clf; plot(x, ux, '-', x, ue);

% shock location
dt = 0.001;
S = 1;
xs = 0;
t = dt:dt:1.0;
xts = zeros(length(t)+1,2);
for i = 1:length(t)
    xm = xs + S*dt;
    tm = i*dt;
    xts(i+1,:) = [xm tm];
    [~, um] = exactsolution(xm, tm, alpha, beta, gamma);        
    S = um/2;
    xs = xm;    
end
figure(1);clf; plot(xts(:,1),xts(:,2));

for i = 1:length(y)     
    [up, um] = exactsolution(x, y(i), alpha, beta, gamma);
    ue(:,i) = 0;
    idx = abs(xts(:,2)-y(i))<1e-8;
    ind = x<=xts(idx,1);
    ue(ind,i) = um(ind);
end
figure(4); clf; plot(x, ux, '-', x, ue);

