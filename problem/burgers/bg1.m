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

UDG01 = initu(mesh,{1.0;0;0});
UH01 = inituhat(master,mesh.elcon,UDG01,1);

S0 = 0.2;
kappa01 = 1.5;
lambda01 = 0.04;

mesh.dgnodes(:,3,:) = lambda01;
[UDG01, UH01] = hdg_solve(master,mesh,app,UDG01,UH01,[]);

figure(1); clf; scaplot(mesh, UDG01(:,1,:),[],2,0); axis on;

eta = 0.8; m = 10;
lambda = ones(m,1)*lambda01;
for i = 2:m
    lambda(i) = lambda(i-1)*eta;
end
kappa = ones(m,1)*kappa01;
for i = 2:m
    kappa(i) = 1 + (kappa(i-1)-1)*eta;
end
[UDG1, UH1, ACG1, mine1] = avloop2(master, mesh, app, UDG01, UH01, S0, lambda, kappa);

% for i=1:length(UDG1)
%     figure(i); clf; scaplot(mesh, UDG1{i}(:,1,:),[],2,0); axis on;
% end
% for i=1:length(UDG1)
%     figure(i); clf; scaplot(mesh, ACG1{i}(:,1,:),[],2,0); axis on;
% end

i = 8; nref=6; alpha = 100;  n = 40; m = 100; polydegree=6; d=2;
div = divergence(UDG1{i}, 1.0);
divmax = max(div(:));
s = limiting(div,0,divmax/2,100,0);
he = meshsize(mesh);
hk =  (kappa(i)*min(he(:)))^2;
s = cgpoisson(mesh, master, s, [hk 1.0]);        
s = s/max(s(:));
a = (s-S0).*(atan(alpha*(s-S0))/pi + 0.5) - atan(alpha)/pi + 0.5;    
S = reshape(a(mesh.t2'), master.npe, 1, mesh.ne);       
figure(1); clf; scaplot(mesh, S(:,1,:),[],2,0); 
[px1, sx1, idx1] = shockpoints(mesh, S, S0, 2*S0, nref);
[pp1, pl1, pr1, pa1] = shockcells(px1, d, n, m);
x = squeeze(pa1(:,1,:))';
y = squeeze(pa1(:,2,:))';
[pn, un] = shocklocation(mesh, S, idx1, x, y, nref+2);
pn(1,:) = [0 0];
pn(end-1:end,:)=xts(end-1:end,:);
poly = polyfit(pn(:,1), pn(:,2), polydegree);
x2 = linspace(0,xts(end,1), 1000);
y2 = polyval(poly,x2);

% figure(1); clf; scaplot(mesh, S(:,1,:),[],2,0); 
% hold on;
% plot(pn(:,1),pn(:,2),'o');
% plot(xts(:,1),xts(:,2) ,'-', 'LineWidth', 2); 
% plot(x1,y1 ,'--', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2); 
% plot(x2,y2,'-.r', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2); 
% 

