function [mesh3,mesh1,mesh2] = mkmesh_unstructuredrect(h,porder)

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

% fit the shock curve with polynomial degree 6
poly = polyfit(pn(:,1) ,pn(:,2), 6);
x1 = linspace(pn(1,1),pn(end,1),8);
y1 = polyval(poly,x1);
xy1 = [x1(:) y1(:)];

% generate the left mesh
pv1 = [-1 0; xy1; -1 1];
[p1,t1]=polymesh({pv1},[1],[0,1],[h,1.3]);

bndexpr = {'all(p(:,2)<1e-3)', ...
           'all(p(:,2)>1-1e-3)',...
           'all(p(:,1)<-1+1e-3)',...
           'true'};                  
mesh1 = mkmesh(p1,t1,porder,bndexpr,0,1);

fb = 4;
mesh1.fcurved = (mesh1.f(:,end)==-fb);
ic = mesh1.fcurved;
mesh1.tcurved = false(size(mesh1.t,1),1);
mesh1.tcurved(mesh1.f(ic,end-1)) = true;
mesh1.dgnodes = makedgnodes(mesh1,@shockdf,fb);

% generate the right mesh
pv2 = [1 0; 1 1; xy1(end:-1:1,:)];
[p2,t2]=polymesh({pv2},[1],[0,1],[h,1.3]);

bndexpr = {'all(p(:,2)<1e-3)', ...
           'all(p(:,2)>1-1e-3)',...
           'all(p(:,1)>1-1e-3)',...
           'true'};                  
mesh2 = mkmesh(p2,t2,porder,bndexpr,0,1);

mesh2.fcurved = (mesh2.f(:,end)==-fb);
ic = mesh2.fcurved;
mesh2.tcurved = false(size(mesh2.t,1),1);
mesh2.tcurved(mesh2.f(ic,end-1)) = true;
mesh2.dgnodes = makedgnodes(mesh2,@shockdf,fb);

% connect the left mesh and the right mesh
[p,t] = connectmesh(p1,t1,p2,t2);
bndexpr = {'all(p(:,2)<min(p0(:,2))+1e-3)','all(p(:,1)>max(p0(:,1))-1e-3)', ...
           'all(p(:,2)>max(p0(:,2))-1e-3)','all(p(:,1)<min(p0(:,1))+1e-3)'};     
mesh3 = mkmesh(p,t,porder,bndexpr,0,1);
mesh3.dgnodes = cat(3,mesh1.dgnodes,mesh2.dgnodes);

% figure(1);clf;meshplot(mesh1,[0 1]);
% hold on;meshplot(mesh2,[0 1]);
% figure(1); clf; simpplot(p,t);
% 
% figure(3); clf; simpplot(p2,t2);
% hold on; simpplot(p1,t1);


