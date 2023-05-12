setapplicationpath('FM/poi');

porder = 5;
ngrid  = 5;
elemtype = 0;
nodetype = 1;
hybrid = 'hdg';

kappa = 1;
c = [0,0]; 
tau = 1;

app.source = 'source';
app.flux = 'flux';
app.fbou = 'fbou';
app.fhat = 'fhat';
app.localsolve=1;
app.arg = {kappa,tau};
app.bcm = [5;5;5;5];
app.bcs = [0;0;0;0];
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

mesh   = mkmesh_square(ngrid,ngrid,porder,0,1,1,elemtype,nodetype);
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);
%mesh = periodic(mesh,{1,'p(:,1)',3,'p(:,1)';4,'p(:,2)',2,'p(:,2)'});

UDG = initu(mesh,{0;0;0});
UH=inituhat(master,mesh.elcon,UDG,1);

% HDG solver
[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,0*UDG);

mesh1 = mkmesh_square(ngrid,ngrid,porder+1,0,1,1,elemtype,nodetype);
master1 = mkmaster(mesh1,2*(porder+1));
[master1,mesh1] = preprocess(master1,mesh1,hybrid);
UDGstar = postprocessnd(master,mesh,master1,mesh1,UDG);

%[UDG,UH]=hdg_solve(master,mesh,app,UDG,UH,SH)
% figure(1); scaplot(mesh,UDG(:,1,:),[],0,1); axis equal; axis tight; colormap jet;
% 
x = (mesh.dgnodes(:,1,:));
y = (mesh.dgnodes(:,2,:));
uex = sin(0.5*pi*x).*sin(0.5*pi*y);
qex = 0*mesh.dgnodes(:,1:2,:);
vex = zeros(master.npv,mesh.nd,mesh.nd,mesh.ne);
qex(:,1,:) = 0.5*pi*cos(0.5*pi*x).*sin(0.5*pi*y);
qex(:,2,:) = 0.5*pi*sin(0.5*pi*x).*cos(0.5*pi*y);
vex(:,1,1,:) = -(0.5*pi)^2*sin(0.5*pi*x).*sin(0.5*pi*y);
vex(:,2,2,:) = -(0.5*pi)^2*sin(0.5*pi*x).*sin(0.5*pi*y);
vex(:,1,2,:) = (0.5*pi)^2*cos(0.5*pi*x).*cos(0.5*pi*y);
vex(:,2,1,:) = (0.5*pi)^2*cos(0.5*pi*x).*cos(0.5*pi*y);

v = UDG(:,1,:);
max(abs(uex(:)-v(:)))
v = -UDG(:,2:3,:);
max(abs(qex(:)-v(:)))

app.kappa = 1.0;
app.tau = tau;
app.coord = 0;
app.param = [kappa,tau];
bcm = app.bcm;
app.fbou = 'ubou';
app.source = 'source';
app.bcm = [5;5;5;5]*0;
%[u,q,uhat,v] = hdgma_poisson(master, mesh, app);
[u,q,uhat,AE,RE,U,W] = hdg_poisson(master, mesh, app);

max(abs(uex(:)-u(:)))
max(abs(qex(:)-q(:)))

% for l=1:2
%   for m=1:2
%     e = vex(:,l,m,:)-v(:,l,m,:);
%     max(abs(e(:)))
%   end
% end
% 
%pause
% 
% figure(1); clf; scaplot(mesh,vex(:,1,1,:),[],0,1); axis equal; axis tight; colormap jet;
% figure(2); clf; scaplot(mesh,v(:,1,1,:),[],0,1); axis equal; axis tight; colormap jet;

[Ae, Re, De, We, Ue, MiB, MiC, MiD, MiE, G, L, H, S] = hdg_precompute(master, mesh, app, UDG);
uh = hdg_linearsystem(Ae, Re, mesh.elcon, 0);
max(abs(uh(:)-uhat(:)))
max(abs(Re(:)-RE(:)))
max(abs(Ae(:)-AE(:)))
max(abs(We(:)-W(:)))
max(abs(Ue(:)-U(:)))

% npv = master.npv;
% npfe = master.npf*size(master.perm,2);
% nd = master.nd;
% ne = size(mesh.dgnodes,3);
% Ba = permute(reshape(Ba,[npv nd npv ne]),[1 3 2 4]);
% Ca = permute(reshape(Ca,[npv nd npfe ne]),[1 3 2 4]);
% max(abs(Ba(:)-B(:)))
max(abs(Ca(:)-C(:)))


return;


nn = [2 4 8 16 32]+1;
for porder = 1:4
    for ii=1:length(nn)
        ngrid = nn(ii);
        poi2d;
        e = calerror(UDG,mesh,master,@exactsol1);          
        erru(ii,porder) = e(1); errq(ii,porder) = sqrt(e(2)^2+e(3)^2);
        errs(ii,porder) = calerror(UDGstar,mesh1,master1,@exactsol1);
    end
end

cu=log(erru(1:end-1,:)./erru(2:end,:))/log(2);
cs=log(errs(1:end-1,:)./errs(2:end,:))/log(2);
cq=log(errq(1:end-1,:)./errq(2:end,:))/log(2);
a=[erru [0 0 0 0; cu]; errq [0 0 0 0; cq]; errs [0 0 0 0; cs]];
a=a(:,[1 5 2 6 3 7 4 8]);
a=[[nn-1 nn-1 nn-1]' a];


figure(1); clf; scaplot3(mesh,UDG(:,1,:),[-0.2 1.2],2); 
%scaplot(mesh,UDG(:,1,:),[0 1],1,1,1)
axis equal; axis tight; colormap jet;
axis([0 1 0 1 -0.2 1.2]);
axis normal; colorbar off;
set(gca,'FontSize',16);
set(gca,'xtick',[0:0.2:1]);
set(gca,'ytick',[0:0.2:1]);
set(gca,'ztick',[-0.2:0.2:1.2]);
xlabel('x','FontSize',18);
ylabel('y','FontSize',18);
box on;

figure(1); clf; scaplot3(mesh1,UDGstar(:,1,:),[-0.2 1.2],2); 
%scaplot(mesh,UDG(:,1,:),[0 1],1,1,1)
axis equal; axis tight; colormap jet;
axis([0 1 0 1 -0.2 1.2]);
axis normal; colorbar off;
set(gca,'FontSize',16);
set(gca,'xtick',[0:0.2:1]);
set(gca,'ytick',[0:0.2:1]);
set(gca,'ztick',[-0.2:0.2:1.2]);
xlabel('x','FontSize',18);
ylabel('y','FontSize',18);
box on;

[RU,RH,RQ] = hdg_residual(master,app,mesh,UDG,UH,0*UDG);

syms x y 
u = sin(0.5*pi*x).*sin(0.5*pi*y);
ux = diff(u,'x');
uxx = diff(ux,'x');
uy = diff(u,'y');
uyy = diff(uy,'y');
f = -simplify(uxx+uyy);

% app.ubou = 'ldgubou';
% UH0 = ldg_uhat(mesh,master,app,UDG(:,1,:));

% app.source = 'source';
% app.flux = 'flux';
% app.fbou = 'ldgfbou';
% app.fhat = 'ldgfhat';
%[RU0,UDG0,UH0] = ldg_residual(master,app,mesh,UDG(:,1,:),0*UDG(:,1,:));

% app.source = 'source';
% app.flux = 'flux';
% app.fbou = 'fbou';
% app.fhat = 'fhat';
% [RU,RH,RQ] = hdg_residual(master,app,mesh,UDG,UH,0*UDG);

% 
app.source = 'source';
app.flux = 'flux';
app.ubou = 'ldgubou';
app.fbou = 'ldgfbou';
app.fhat = 'ldgfhat';
tol = 1e-8;
dt = 1e-1*ones(1000,1);
[UDGA,UHA,normR] = ardm(master,app,mesh,0*UDG(:,1,:),dt,tol);

figure(1); scaplot(mesh,UDGA(:,1,:),[],0,1); axis equal; axis tight; colormap jet;

% [un,normR] = restartedardm(master,app,mesh,UDG,UH,0*SH,dt,tol);


% % HDG postprocessing 
% mesh1 = mkmesh_square(ngrid,ngrid,porder+1,0,1,1,elemtype,nodetype);
% master1 = mkmaster(mesh1,2*(porder+1));
% [master1,mesh1] = preprocess(master1,mesh1,hybrid);
% UDGstar = postprocessnd(master,mesh,master1,mesh1,UDG);
% 
% figure(2); scaplot(mesh1,UDGstar(:,1,:),[],0,1); axis equal; axis tight;
% %figure(2); scaplot(mesh,u(:,1,:),[],0,1); axis equal; axis tight;
% 
% x = (mesh1.dgnodes(:,1,:));
% y = (mesh1.dgnodes(:,2,:));
% u = sin(pi*x).*sin(pi*y);
% v = UDGstar(:,1,:);
% max(abs(u(:)-v(:)))
% 
