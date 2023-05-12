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

% Uniform grid
mesh = mkmesh_rect(ngrid*2+1,ngrid+1,porder,0,[-1 1 0 1],1,1);
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);
mesh = mkcgmesh(mesh);
mesh.ib = [];
mesh.in = 1:size(mesh.p2,1);

UDG0 = initu(mesh,{1.0;0;0});
av = [0.02 0.4 0.3 0.25 0.2 0.15 0.1];
hk = 0.005;
href = 1/(ngrid*porder);
alpha = 100;
[UDG, UH, ACG, mine] = avloop(master, mesh, app, UDG0, av, href, hk, alpha);

for i = 1:length(av)
    figure(i); clf;scaplot(mesh,UDG{i}(:,1,:),[],2,0); axis equal; axis tight;
end

[~,ind1] = min(abs(mine-0.5));
UDG1 = UDG{ind1};
UH1 = UH{ind1};
ACG1 = ACG{ind1-1};
figure(1); clf;scaplot(mesh,UDG1(:,1,:),[0 2],2,0); axis equal; axis tight;
figure(2); clf;scaplot(mesh,ACG1(:,1,:),[],2,0); axis equal; axis tight;
figure(3); clf;scaplot(mesh,UDG1(:,1,:),[0 2],2,1); axis equal; axis tight;
figure(4); clf;scaplot(mesh,ACG1(:,1,:),[],2,1); axis equal; axis tight;

[cgnodes,cgelcon,rowent2elem,colent2elem,cgent2dgent] = mkcgent2dgent(mesh.dgnodes(:,1:2,:),1e-8);
s = ucg2udg(udg2ucg(divergence(UDG1,href), cgent2dgent, rowent2elem), cgelcon);
figure(5); clf;scaplot(mesh,s,[],2,1); axis equal; axis tight;

s0 = 0.1; s1 = 0.18; nref = 6; d = 2; n = 20; m = 1000; polydegree=6; 
[px, sx, idx] = shockpoints(mesh, s, s0, s1, nref);
[pp, pl, pr, pa] = shockcells(px, d, n, m);
x = squeeze(pa(:,1,:))';
y = squeeze(pa(:,2,:))';
[pn, un] = shocklocation(mesh, s, idx, x, y, nref);

%mesh2 = mkmesh_shock2(porder);
mesh2 = mkmesh_shock1(porder, pn, polydegree, 18);
master = mkmaster(mesh2,2*porder);
[master,mesh2] = preprocess(master,mesh2,hybrid);
mesh2 = mkcgmesh(mesh2);
mesh2.ib = [];
mesh2.in = 1:size(mesh2.p2,1);

UDG0 = initu(mesh2,{0.0;0;0});
av = [0.01 [1 1/2 1/4 1/8 1/16 1/32 1/64 1/128 1/256 1/1024]*0.0025];
hk = 0.000001;
[UDG, UH, ACG, mine] = avloop(master, mesh2, app, UDG0, av, href, hk, alpha);

for i = 1:length(av)
    figure(i); clf;scaplot(mesh2,UDG{i}(:,1,:),[],2,0); axis equal; axis tight;
end
for i = 1:length(av)
    figure(i); clf;scaplot(mesh2,ACG{i}(:,1,:),[],2,0); axis equal; axis tight;
end

ind2 = 10;
UDG2 = UDG{ind2};
UH2 = UH{ind2};
ACG2 = ACG{ind2-1};

x = linspace(-1, 1, 2000);
y = [0.05 0.1:0.1:0.9];
ux1 = zeros(length(x), length(y));
for i = 1:length(y)
    xy = [x(:) y(i)*ones(length(x),1)];
    ux1(:,i) = fieldatx(mesh,UDG1(:,1,:),xy,20);
end
figure(4); clf; plot(x, ux, '-');

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
    [~, um] = exactsolution(xm, tm, 2, 4, 2);        
    S = um/2;
    xs = xm;    
end
figure(1);clf; plot(xts(:,1),xts(:,2));

ue = zeros(length(x), length(y));
for i = 1:length(y)     
    [up, um] = exactsolution(x, y(i), 2, 4, 2);
    ue(:,i) = 0;
    idx = abs(xts(:,2)-y(i))<1e-8;
    ind = x<=xts(idx,1);
    ue(ind,i) = um(ind);
end
figure(4); clf; plot(x, ux, '-', x, ue);

