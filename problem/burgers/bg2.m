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
mesh2 = mkmesh_shock1(porder, pn, polydegree, 18);
master = mkmaster(mesh2,2*porder);
[master,mesh2] = preprocess(master,mesh2,hybrid);
mesh2 = mkcgmesh(mesh2);
mesh2.ib = [];
mesh2.in = 1:size(mesh2.p2,1);

UDG02 = initu(mesh2,{0.0;0;0});
UH02 = inituhat(master,mesh2.elcon,UDG02,1);

S0 = 0.2;
kappa02 = 5;
lambda02 = 0.01;

mesh2.dgnodes(:,3,:) = lambda02;
[UDG02, UH02] = hdg_solve(master,mesh2,app,UDG02,UH02,[]);
figure(1); clf; scaplot(mesh2, UDG02(:,1,:),[],2,0); axis on;

eta = 0.7; m = 12;
lambda2 = ones(m,1)*lambda02;
for i = 2:m
    lambda2(i) = lambda2(i-1)*eta;
end
kappa2 = ones(m,1)*kappa02;
for i = 2:m
    kappa2(i) = 1 + (kappa2(i-1)-1)*eta;
end
[UDG2, UH2, ACG2, mine2] = avloop2(master, mesh2, app, UDG02, UH02, S0, lambda2, kappa2);

% 
for i=1:length(UDG2)
    figure(i); clf; scaplot(mesh2, UDG2{i}(:,1,:),[],2,0); axis on;
end
for i=1:length(UDG2)
    figure(i); clf; scaplot(mesh2, ACG2{i}(:,1,:),[],2,0); axis on;
end

x = linspace(-1, 1, 2000);
y = [0.05 0.1:0.1:0.9];
ux3 = zeros(length(x), length(y));
for i = 1:length(y)
    i
    xy = [x(:) y(i)*ones(length(x),1)];
    ux3(:,i) = fieldatx(mesh,UDG1{4}(:,1,:),xy,16);
end

x = linspace(-1, 1, 2000);
y = [0.05 0.1:0.1:0.9];
ux4 = zeros(length(x), length(y));
for i = 1:length(y)
    i
    xy = [x(:) y(i)*ones(length(x),1)];
    ux4(:,i) = fieldatx(mesh2,UDG2{12}(:,1,:),xy,16);
end

ii = [1 3 5 10];
figure(6); clf; hold on;
for i = 1:length(ii)
    plot(x, ue(:,ii(i)), '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
    plot(x, ux3(:,ii(i)), '-.', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2); 
    plot(x, ux4(:,ii(i)), '--', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2); 
end
set(gca,'FontSize',18); 
set(gca,'LooseInset',get(gca,'TightInset'))
text(0.065,1.825,'$\leftarrow t = 0.05$', 'interpreter', 'latex', 'FontSize',18);
text(0.185,1.505,'$\leftarrow t = 0.20$', 'interpreter', 'latex', 'FontSize',18);
text(0.325,1.275,'$\leftarrow t = 0.40$', 'interpreter', 'latex', 'FontSize',18);
text(0.605,0.975,'$\leftarrow t = 0.90$', 'interpreter', 'latex', 'FontSize',18);
xlabel("$x$", 'interpreter', 'latex', 'FontSize', 22);
ylabel("$u$", 'interpreter', 'latex', 'FontSize', 22);
box on;
leg = legend({'Exact solution', 'Regular mesh', 'Shock-aligned mesh'}, 'FontSize', 16, 'Location', 'NW');
leg.ItemTokenSize = [40,10];

x = linspace(-1, 1, 2000);
y = [0.05 0.1:0.1:0.9];
ux = zeros(length(x), length(y), 12);
for m = 1:12
  for i = 1:length(y)
      xy = [x(:) y(i)*ones(length(x),1)];
      ux(:,i,m) = fieldatx(mesh2,UDG2{m}(:,1,:),xy,18);
  end
end

figure(6); clf; hold on;
set(gca,'FontSize',18); 
plot(x, ue(:,i), '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 3); 
plot(x, squeeze(ux(:,i,[1 2 4 8 10])), '-.', 'LineWidth', 3); 
xlabel("$x$", 'interpreter', 'latex', 'FontSize', 24);
ylabel("$u$", 'interpreter', 'latex', 'FontSize', 24);
box on;
leg = legend({'Exact solution', 'n=1', 'n=2', 'n=4', 'n=8', 'n=12'}, 'FontSize', 18, 'Location', 'NW');
leg.ItemTokenSize = [30,10];
ax = gca;
fn = "bg_pro" + num2str(i) + ".png";
exportgraphics(ax,fn,'Resolution',200); 


set(gca,'FontSize',22); 
xlabel("$x$", 'interpreter', 'latex', 'FontSize', 28);
ylabel("$u$", 'interpreter', 'latex', 'FontSize', 28);
legend('off');


