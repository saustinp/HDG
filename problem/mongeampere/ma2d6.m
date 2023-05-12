setapplicationpath('FM/poi');

porder = 3;
ngrid  = 81;
elemtype = 1;
nodetype = 1;
hybrid = 'hdg';

app.tau = 1;
app.bcm = [1;3;1;3];
app.bcd = [1 1; 1 1; 1 1; 1 1];
app.fbou = 'ubouma7';
app.source = 'sourcema6';

a = 1;
b = 2;
c = 4;
t1 = pi/2;
t2 = 3*pi/2;
mesh = mkmesh_square(ngrid,ngrid,porder,1,1,1,1,1);
mesh = mkmesh_halfcircleellipse(mesh,a,b,c,t1,t2);
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);
 
x = (mesh.dgnodes(:,1,:));
y = (mesh.dgnodes(:,2,:));
r = sqrt((x-1.5).^2+y.^2);
a1 = 10;
a2 = 6;
a = 3;
rho = 1 + a1*sech(a2*(r.^2-a^2));
L = averagevector(master,mesh);
theta = sum(L(:).*rho(:))/sum(L(:));

figure(1); clf; scaplot(mesh,rho,[],0,1); axis equal; axis tight; colormap default;
figure(2); clf; meshplot(mesh); axis equal;

u = (x.^2+y.^2)/2;
uavg = sum(L(:).*u(:));
u = u - uavg;
app.param = [theta a1 a2 a uavg];
uhat = inituhat(master,mesh.elcon,u,1);
uhat = uhat(:);

app.neumman = 1;
[u,q,uhat,v] = hdg_ma5(master, mesh, app, u, uhat);

t=linspace(0,2*pi,1000);
x=cos(t);
y=sin(t);

figure(2); clf; meshplot(mesh,1);   
hold on;
plot(x,y,'-r');

mesh2 = mesh; 
mesh2.dgnodes = q;
figure(3); clf; meshplot(mesh2,1);
hold on;
plot(x,y,'-r');

%checkneumann(master, mesh, app, q);

%figure(1); clf; scaplot(mesh,v(:,1,1,:),[],0,1); axis equal; axis tight; colormap jet;


% x = (mesh.dgnodes(:,1,:));
% y = (mesh.dgnodes(:,2,:));
% r2 = x.^2 + y.^2;
% rho = (r2 - 1).^3  - x.^2.*y.^3;
% rho = 1 + exp(-10*rho.^2);
% figure(1); clf; meshplot(mesh,1); axis equal;
% 
% 
% mesh   = mkmesh_rect(ngrid,ngrid,porder,0,4*[-0.5 0.5 -0.5 0.5],elemtype,nodetype);
% L = averagevector(master,mesh);
% theta = sum(L(:).*rho(:))/sum(L(:));
% 
% u = (x.^2+y.^2)/2;
% uavg = sum(L(:).*u(:));
% u = u - uavg;
% app.param = [theta a1 a2 a uavg];
% uhat = inituhat(master,mesh.elcon,u,1);
% uhat = uhat(:);
% 
% figure(1); clf; scaplot(mesh,rho,[],0,1); axis equal; axis tight; colormap jet;
% %figure(2); clf; scaplot(mesh,u,[],0,1); axis equal; axis tight; colormap jet;
% 
% app.neumman = 1;
% [u,q,uhat,v] = hdg_ma4(master, mesh, app, u, uhat);
% 
% figure(1); clf; scaplot(mesh,v(:,1,1,:),[],0,1); axis equal; axis tight; colormap jet;

return;

