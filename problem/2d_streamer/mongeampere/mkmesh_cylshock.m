function mesh = mkmesh_cylshock(porder, pn)

m = 31;
n = 10;
a = 1; % inner radius
t1 = pi/2;
t2 = 3*pi/2;

tm = t2 + linspace(0,1,m)'*(t1-t2);
ds = [0 0.001 0.003 0.007 0.016 0.032 0.064 0.128];
ds = [-ds(end:-1:2) ds];
ds = [ds 0.21 0.4 0.7 1.1734];

rn = sqrt(pn(:,1).^2 + pn(:,2).^2);
tn = atan(pn(:,2)./pn(:,1)) + pi;
polyn = polyfit(tn, rn, 6);

[pa,ta,dga] = cylshockgrid(tm, ds, polyn, porder);

nx = cos(tm);
ny = sin(tm);
dm = polyval(polyn, tm); 
xm = dm.*cos(tm);
ym = dm.*sin(tm);
xn = xm + ds(1)*nx;
yn = ym + ds(1)*ny;

[p,t] = squaremesh(n,m,0,1);
mesh = mkmesh(p,t,porder,{'true'},1,1);
%dgb = mesh.dgnodes;

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

[p, t] = connectmesh(pa, ta, mesh.p, mesh.t, 1e-5);
dgnodes = cat(3, dga, mesh.dgnodes);

bndexpr = {'all(p(:,1)>-1e-6)','all(sqrt(sum(p.^2,2))<2)','all(sqrt(sum(p.^2,2))>2)'};  
mesh = mkmesh(p,t,porder,bndexpr,1,1);
mesh.dgnodes = dgnodes;

figure(1); clf; 
meshplot(mesh,[1 0]);
axis equal; axis tight; axis on;

figure(2); clf; hold on;
boundaryplot(mesh,0,1);
boundaryplot(mesh,0,2);
boundaryplot(mesh,0,3);
axis equal; axis tight; axis on;
