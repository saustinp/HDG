function mesh = mkmesh_shockrect(porder)

pn = [0                   0
   0.040000000000000   0.050000000000000
   0.085000000000000   0.100000000000000
   0.125000000000000   0.150000000000000
   0.165000000000000   0.200000000000000
   0.200000000000000   0.250000000000000
   0.237500000000000   0.300000000000000
   0.265000000000000   0.350000000000000
   0.300000000000000   0.400000000000000
   0.335000000000000   0.450000000000000
   0.362500000000000   0.500000000000000
   0.390000000000000   0.550000000000000
   0.417500000000000   0.600000000000000
   0.445000000000000   0.650000000000000
   0.475000000000000   0.700000000000000
   0.500000000000000   0.750000000000000
   0.525000000000000   0.800000000000000
   0.552500000000000   0.850000000000000
   0.577500000000000   0.900000000000000
   0.600000000000000   0.950000000000000
   0.625000000000000   1.000000000000000];

pfit = 6; n = 12;

% fit the shock curve with polynomial degree 6
poly = polyfit(pn(:,1) ,pn(:,2), pfit);
x1 = linspace(pn(1,1),pn(end,1),n);
y1 = polyval(poly,x1);
y1 = max(y1,0);
y1 = min(y1,1);
y1(1) = 0;
%xy1 = [x1(:) y1(:)];

ds = [0 0.0015 0.003 0.006 0.012 0.024 0.048 0.096 0.155];
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
           'all((p(:,1)+0.1)/0.68 < p(:,2))',...
           'true'};                  
meshl = mkmesh(pl,tl,porder,bndexpr,1,1);
% fb = 4;
% meshl.fcurved = (meshl.f(:,end)==-fb);
% ic = meshl.fcurved;
% meshl.tcurved = false(size(meshl.t,1),1);
% meshl.tcurved(meshl.f(ic,end-1)) = true;
% meshl.dgnodes = makedgnodes(meshl,@shockdf,fb);

bndexpr = {'all(p(:,2)<1e-3)', ...
           'all(p(:,2)>1-1e-3)',...
           'all((p(:,1)-0.1)/0.6 < p(:,2))',...
           'true'};                         
meshr = mkmesh(pr,tr,porder,bndexpr,1,1);
% fb = 3;
% meshr.fcurved = (meshr.f(:,end)==-fb);
% ic = meshr.fcurved;
% meshr.tcurved = false(size(meshr.t,1),1);
% meshr.tcurved(meshr.f(ic,end-1)) = true;
% meshr.dgnodes = makedgnodes(meshr,@shockdf,fb);

pbl = boundarypoints(meshl.p,meshl.f,3);
pbr = boundarypoints(meshr.p,meshr.f,4);
figure(3);clf;plot(pbl(:,1),pbl(:,2),'-o',pbr(:,1),pbr(:,2),'-o');

k = [8 10 7];
x2 = linspace(pbl(end,1), -1.0, k(1))';
y2 = ones(k(1),1);
x3 = -ones(k(2),1);
y3 = linspace(1.0, 0, k(2))';
x4 = linspace(-1.0, pbl(1,1), k(3))';
y4 = zeros(k(3),1);
pv = [x2(2:end) y2(2:end); x3(2:end) y3(2:end); x4(2:end-1) y4(2:end-1)];
pbl = [pbl; pv];

k = [3 10 9];
% x2 = linspace(pbr(end,1),1.0,k(1))';
% y2 = ones(k(1),1);
% x3 = ones(k(2),1);
% y3 = linspace(1.0, 0, k(2))';
% x4 = linspace(1.0, pbr(1,1), k(3))';
% y4 = zeros(k(3),1);
% pv = [x2(2:end) y2(2:end); x3(2:end) y3(2:end); x4(2:end-1) y4(2:end-1)];
% pbr = [pbr; pv];
x2 = linspace(pbr(1,1),1.0, k(3))';
y2 = zeros(k(3),1);
x3 = ones(k(2),1);
y3 = linspace(0, 1.0, k(2))';
x4 = linspace(1.0,pbr(end,1),k(1))';
y4 = ones(k(1),1);
pv = [x2(2:end) y2(2:end); x3(2:end) y3(2:end); x4(2:end-1) y4(2:end-1)];
pbr = [pv; pbr(end:-1:1,:)];

poly2gmsh('leftdomain.geo', pbl, 0.3);
gmshmatlab('leftdomain', '-2');
[p1,t1] = gmsh2pq('leftdomain.msh');

poly2gmsh('rightdomain.geo', pbr, 0.3);
gmshmatlab('rightdomain', '-2');
[p2,t2] = gmsh2pq('rightdomain.msh');

figure(2);clf;simpplot(p1,t1); hold on; plot(pbl(:,1),pbl(:,2),'-o'); axis on;
figure(3);clf;simpplot(p2,t2); hold on; plot(pbr(:,1),pbr(:,2),'-o'); axis on;

[p,t] = connectmesh(p1,t1,pl,tl,1e-4);
[p,t] = connectmesh(p,t,pr,tr,1e-4);
[p,t] = connectmesh(p,t,p2,t2,1e-4);
bndexpr = {'all(p(:,2)<min(p0(:,2))+1e-3)','all(p(:,1)>max(p0(:,1))-1e-3)', ...
           'all(p(:,2)>max(p0(:,2))-1e-3)','all(p(:,1)<min(p0(:,1))+1e-3)'};     
mesh = mkmesh(p,t,porder,bndexpr,1,1);
figure(1);clf; meshplot(mesh,[1 0]); axis on;



% figure(2);clf;meshplot(meshl,[1 1]);
% hold on;meshplot(meshr,[1 1]);
% axis on;
% 
% % generate the left mesh
% pv1 = [-1 0; xy1; -1 1];
% [p1,t1]=polymesh({pv1},[1],[0,1],[h,1.3]);
% 
% bndexpr = {'all(p(:,2)<1e-3)', ...
%            'all(p(:,2)>1-1e-3)',...
%            'all(p(:,1)<-1+1e-3)',...
%            'true'};                  
% mesh1 = mkmesh(p1,t1,porder,bndexpr,0,1);
% 
% fb = 4;
% mesh1.fcurved = (mesh1.f(:,end)==-fb);
% ic = mesh1.fcurved;
% mesh1.tcurved = false(size(mesh1.t,1),1);
% mesh1.tcurved(mesh1.f(ic,end-1)) = true;
% mesh1.dgnodes = makedgnodes(mesh1,@shockdf,fb);
% 
% % generate the right mesh
% pv2 = [1 0; 1 1; xy1(end:-1:1,:)];
% [p2,t2]=polymesh({pv2},[1],[0,1],[h,1.3]);
% 
% bndexpr = {'all(p(:,2)<1e-3)', ...
%            'all(p(:,2)>1-1e-3)',...
%            'all(p(:,1)>1-1e-3)',...
%            'true'};                  
% mesh2 = mkmesh(p2,t2,porder,bndexpr,0,1);
% 
% mesh2.fcurved = (mesh2.f(:,end)==-fb);
% ic = mesh2.fcurved;
% mesh2.tcurved = false(size(mesh2.t,1),1);
% mesh2.tcurved(mesh2.f(ic,end-1)) = true;
% mesh2.dgnodes = makedgnodes(mesh2,@shockdf,fb);
% 
% % connect the left mesh and the right mesh
% [p,t] = connectmesh(p1,t1,p2,t2);
% bndexpr = {'all(p(:,2)<min(p0(:,2))+1e-3)','all(p(:,1)>max(p0(:,1))-1e-3)', ...
%            'all(p(:,2)>max(p0(:,2))-1e-3)','all(p(:,1)<min(p0(:,1))+1e-3)'};     
% mesh3 = mkmesh(p,t,porder,bndexpr,0,1);
% mesh3.dgnodes = cat(3,mesh1.dgnodes,mesh2.dgnodes);
% 
% figure(1);clf;meshplot(mesh1,[1 0]);
% hold on;meshplot(mesh2,[1 0]);
% %figure(1); clf; simpplot(p,t);
% % 
% % figure(3); clf; simpplot(p2,t2);
% % hold on; simpplot(p1,t1);
% 
% 
