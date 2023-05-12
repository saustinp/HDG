function mesh = mkmesh_shock1(porder, pn, pfit, n)

% fit the shock curve with polynomial degree 6
poly = polyfit(pn(:,1) ,pn(:,2), pfit);
x1 = linspace(pn(1,1),pn(end,1),n);
y1 = polyval(poly,x1);
y1 = max(y1,0);
y1 = min(y1,1);
y1(1) = 0;
%xy1 = [x1(:) y1(:)]

ds = [0 0.0008 0.002 0.004 0.008 0.016 0.032 0.064 0.128];
m = length(ds);
sx = zeros(n, m);
sy = 0*sx;
for i = 1:m
    sx(:,i) = x1(:) + ds(i);
    sy(:,i) = y1(:);
end
sx=sx';
sy=sy';
pr = [sx(:) sy(:)];
t = [1 2 m+2 m+1];
t = kron(t,ones(n-1,1))+kron(ones(size(t)),(0:n-2)'*m);
tr = kron(t,ones(m-1,1))+kron(ones(size(t)),(0:m-2)');

%ds = [0 0.005 0.015 0.035 0.075 0.155];
ds = -ds(end:-1:1);
m = length(ds);
sx = zeros(n, m);
sy = 0*sx;
for i = 1:m
    sx(:,i) = x1(:) + ds(i);
    sy(:,i) = y1(:);
end
sx=sx';
sy=sy';
pl = [sx(:) sy(:)];
t = [1 2 m+2 m+1];
t = kron(t,ones(n-1,1))+kron(ones(size(t)),(0:n-2)'*m);
tl = kron(t,ones(m-1,1))+kron(ones(size(t)),(0:m-2)');

bndexpr = {'all(p(:,2)<1e-3)', ...
           'all(p(:,2)>1-1e-3)',...
           'all((p(:,1)+2e-4)/0.57 < p(:,2))',...
           'true'};                  
meshl = mkmesh(pl,tl,porder,bndexpr,1,1);

bndexpr = {'all(p(:,2)<1e-3)', ...
           'all(p(:,2)>1-1e-3)',...
           'all((p(:,1)-0.1)/0.6 < p(:,2))',...
           'true'};                         
meshr = mkmesh(pr,tr,porder,bndexpr,1,1);

pbl = boundarypoints(meshl.p,meshl.f,3);
pbr = boundarypoints(meshr.p,meshr.f,4);

k = [14 10 10];
x2 = logdec(linspace(pbl(end,1), -1.0, k(1)),0.5)';
y2 = ones(k(1),1);
x3 = -ones(k(2),1);
y3 = linspace(1.0, 0, k(2))';
x4 = logdec(linspace(-1.0, pbl(1,1), k(3))',0.9);
y4 = zeros(k(3),1);
pv = [x2(2:end) y2(2:end); x3(2:end) y3(2:end); x4(2:end-1) y4(2:end-1)];
pbl = [pbl; pv];

%k = [4 12 12];
k = [4 10 10];
x2 = linspace(pbr(end,1),1.0,k(1))';
y2 = ones(k(1),1);
x3 = ones(k(2),1);
y3 = linspace(1.0, 0, k(2))';
x4 = loginc(linspace(1.0, pbr(1,1), k(3)),1.0)';
y4 = zeros(k(3),1);
pv = [x2(2:end) y2(2:end); x3(2:end) y3(2:end); x4(2:end-1) y4(2:end-1)];
pbr = [pbr; pv];

poly2gmsh('leftdomain.geo', pbl, 0.25);
gmshmatlab('leftdomain', '-2');
[p1,t1] = gmsh2pq('leftdomain.msh');

poly2gmsh('rightdomain.geo', pbr, 0.25);
gmshmatlab('rightdomain', '-2');
[p2,t2] = gmsh2pq('rightdomain.msh');

% figure(1);clf; meshplot(meshl); hold on; simpplot(p1,t1); axis on;
% figure(2);clf; meshplot(meshr); hold on; simpplot(p2,t2); axis on;
figure(1);clf;simpplot(p1,t1); hold on; plot(pbl(:,1),pbl(:,2),'-o');
figure(2);clf;simpplot(p2,t2); hold on; plot(pbr(:,1),pbr(:,2),'-o');

[p,t] = connectmesh(p1,t1,pl,tl,1e-4);
[p,t] = connectmesh(p,t,pr,tr,1e-4);
[p,t] = connectmesh(p,t,p2,t2,1e-4);
bndexpr = {'all(p(:,2)<min(p0(:,2))+1e-3)','all(p(:,1)>max(p0(:,1))-1e-3)', ...
           'all(p(:,2)>max(p0(:,2))-1e-3)','all(p(:,1)<min(p0(:,1))+1e-3)'};     
mesh = mkmesh(p,t,porder,bndexpr,1,1);

figure(3);clf; meshplot(mesh,[1 0]);

