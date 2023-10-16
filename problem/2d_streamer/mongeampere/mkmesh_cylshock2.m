function mesh = mkmesh_cylshock2(porder, pn)

m = 31;
n = 14;
k = 8;
a = 1; % inner radius
t1 = pi/2;
t2 = 3*pi/2;

tm = t2 + linspace(0,1,m)'*(t1-t2);
ds = [0 0.001 0.003 0.007 0.016 0.032 0.064 0.128  0.21 0.4 0.7 1.1734];

rn = sqrt(pn(:,1).^2 + pn(:,2).^2);
tn = atan(pn(:,2)./pn(:,1)) + pi;
polyn = polyfit(tn, rn, 6);

% inner mesh
nx = cos(tm);
ny = sin(tm);
dm = polyval(polyn, tm); 
xm = dm.*cos(tm);
ym = dm.*sin(tm);
xn = xm + ds(1)*nx;
yn = ym + ds(1)*ny;

[p,t] = squaremesh(n,m,0,1);
p(:,1) = loginc(logdec(p(:,1),2),4);
%p(:,1) = loginc(p(:,1),1.5);
mesh = mkmesh(p,t,porder,{'true'},1,1);
mesh.dgnodes(:,1,:) = loginc(logdec(mesh.dgnodes(:,1,:),2),4);
%mesh.dgnodes(:,1,:) = loginc(mesh.dgnodes(:,1,:),1.5);

polym = polyfit(tm, sqrt(xn.^2 + yn.^2), 6);
tb = t2 + p(:,2)*(t1-t2);
d = polyval(polym, tb); % outer radius
r = d + p(:,1).*(a-d);
p(:,1) = r.*cos(tb);
p(:,2) = r.*sin(tb);

tb = t2 + mesh.dgnodes(:,2,:)*(t1-t2);
d = polyval(polym, tb); % outer radius
r = d + mesh.dgnodes(:,1,:).*(a-d);
mesh.dgnodes(:,1,:) = r.*cos(tb);
mesh.dgnodes(:,2,:) = r.*sin(tb);
mesh.p = p;
mesh1 = mesh;

figure(1); clf; 
meshplot(mesh,[1 0]);
axis equal; axis tight; axis on;

% figure(2); clf; hold on;
% boundaryplot(mesh,0,1);
% boundaryplot(mesh,0,2);
% boundaryplot(mesh,0,3);
% axis equal; axis tight; axis on;

% outer mesh
xk = xm + ds(end)*nx;
yk = ym + ds(end)*ny;

[p,t] = squaremesh(k,m,0,1);
p(:,1) = logdec(p(:,1),3.5);
mesh = mkmesh(p,t,porder,{'true'},1,1);
mesh.dgnodes(:,1,:) = logdec(mesh.dgnodes(:,1,:),3.5);

polyk = polyfit(tm, sqrt(xk.^2 + yk.^2), 6);
tb = t2 + p(:,2)*(t1-t2);
d = polyval(polym, tb); % outer radius
dk = polyval(polyk, tb); % outer radius
r = dk + p(:,1).*(d-dk);
p(:,1) = r.*cos(tb);
p(:,2) = r.*sin(tb);

tb = t2 + mesh.dgnodes(:,2,:)*(t1-t2);
d = polyval(polym, tb); % outer radius
dk = polyval(polyk, tb); % outer radius
r = dk + mesh.dgnodes(:,1,:).*(d-dk);
mesh.dgnodes(:,1,:) = r.*cos(tb);
mesh.dgnodes(:,2,:) = r.*sin(tb);
mesh.p = p;
mesh2 = mesh;

figure(2); clf; 
meshplot(mesh,[1 0]);
axis equal; axis tight; axis on;

[p, t] = connectmesh(mesh1.p, mesh1.t, mesh2.p, mesh2.t, 1e-5);
dgnodes = cat(3, mesh1.dgnodes, mesh2.dgnodes);

bndexpr = {'all(p(:,1)>-1e-6)','all(sqrt(sum(p.^2,2))<2)','all(sqrt(sum(p.^2,2))>2)'};  
mesh = mkmesh(p,t,porder,bndexpr,1,1);
mesh.dgnodes = dgnodes;

m = 1000;
tm = t2 + linspace(0,1,m)'*(t1-t2);
dm = polyval(polyn, tm); 
xm = dm.*cos(tm);
ym = dm.*sin(tm);

figure(3); clf; 
meshplot(mesh,[1 0]);
hold on;
%plot(pn(:,1),pn(:,2),'-r','LineWidth',3);
%plot(xm,ym,'-r','LineWidth',3);
axis equal; axis tight; axis on;

figure(4); clf; hold on;
boundaryplot(mesh,0,1);
boundaryplot(mesh,0,2);
boundaryplot(mesh,0,3);

[mesh1.ne mesh2.ne mesh.ne]