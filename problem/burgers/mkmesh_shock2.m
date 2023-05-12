function mesh = mkmesh_shock2(porder)

pn = [0.000001215501216  -0.000000000000000
   0.048484591734592   0.050000000000000
   0.091305682305682   0.100000000000000
   0.132035282035282   0.150000000000000
   0.172230158730159   0.200000000000000
   0.209386464386464   0.250000000000000
   0.243709030459030   0.300000000000000
   0.276974402974403   0.350000000000000
   0.309913507663508   0.400000000000000
   0.340947944697945   0.450000000000000
   0.372069751569751   0.500000000000000
   0.400867419367419   0.550000000000000
   0.429746535496535   0.600000000000000
   0.458430661180661   0.650000000000000
   0.484855316355316   0.700000000000000
   0.511134117884118   0.750000000000000
   0.537503773253773   0.800000000000000
   0.563660660660661   0.850000000000000
   0.587777368277368   0.900000000000000
   0.612335744835745   0.950000000000000
   0.636274768274768   1.000000000000000];

pfit = 6; n = 18;

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
% figure(1);clf;meshplot(meshl);
% hold on; plot([2e-4 0.575], [0 1]);

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
%figure(3);clf;plot(pbl(:,1),pbl(:,2),'-o',pbr(:,1),pbr(:,2),'-o');

k = [18 14 12];
x2 = linspace(pbl(end,1), -1.0, k(1))';
y2 = ones(k(1),1);
x3 = -ones(k(2),1);
y3 = linspace(1.0, 0, k(2))';
x4 = linspace(-1.0, pbl(1,1), k(3))';
y4 = zeros(k(3),1);
pv = [x2(2:end) y2(2:end); x3(2:end) y3(2:end); x4(2:end-1) y4(2:end-1)];
pbl = [pbl; pv];

k = [4 12 12];
x2 = linspace(pbr(end,1),1.0,k(1))';
y2 = ones(k(1),1);
x3 = ones(k(2),1);
y3 = linspace(1.0, 0, k(2))';
x4 = linspace(1.0, pbr(1,1), k(3))';
y4 = zeros(k(3),1);
pv = [x2(2:end) y2(2:end); x3(2:end) y3(2:end); x4(2:end-1) y4(2:end-1)];
pbr = [pbr; pv];

poly2gmsh('leftdomain.geo', pbl, 0.22);
gmshmatlab('leftdomain', '-2');
[p1,t1] = gmsh2pq('leftdomain.msh');

poly2gmsh('rightdomain.geo', pbr, 0.25);
gmshmatlab('rightdomain', '-2');
[p2,t2] = gmsh2pq('rightdomain.msh');

figure(1);clf; meshplot(meshl); hold on; simpplot(p1,t1); axis on;
figure(2);clf; meshplot(meshr); hold on; simpplot(p2,t2); axis on;

[p,t] = connectmesh(p1,t1,pl,tl,1e-4);
[p,t] = connectmesh(p,t,pr,tr,1e-4);
[p,t] = connectmesh(p,t,p2,t2,1e-4);
bndexpr = {'all(p(:,2)<min(p0(:,2))+1e-3)','all(p(:,1)>max(p0(:,1))-1e-3)', ...
           'all(p(:,2)>max(p0(:,2))-1e-3)','all(p(:,1)<min(p0(:,1))+1e-3)'};     
mesh = mkmesh(p,t,porder,bndexpr,1,1);
figure(1);clf; meshplot(mesh,[1 0]);

figure(2);clf;simpplot(p1,t1); hold on; plot(pbl(:,1),pbl(:,2),'-o');
figure(3);clf;simpplot(p2,t2); hold on; plot(pbr(:,1),pbr(:,2),'-o');
