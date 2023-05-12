S = [0 0; 1 0; 0 1; 1 1];
L = [-1 0; 0 0; -1 1; 0.5 1];
R = [0 0; 1.0 0; 0.5 1; 1 1];
xrect = [-1 1 0 1];

m = 20; n = m/2; 
[p,t] = squaremesh(2*m+1,n+1,0,0);
p(:,1) = xrect(1) + (xrect(2)-xrect(1))*p(:,1);
p(:,2) = xrect(3) + (xrect(4)-xrect(3))*p(:,2);

bndexpr = {'all(p(:,2)<min(p0(:,2))+1e-3)','all(p(:,1)>max(p0(:,1))-1e-3)', ...
           'all(p(:,2)>max(p0(:,2))-1e-3)','all(p(:,1)<min(p0(:,1))+1e-3)'};     
mesh = mkmesh(p,t,porder,bndexpr,0,1);

% element centers
nd = mesh.nd;
[ne,nv] = size(mesh.t);
p = reshape(mesh.p(mesh.t',:),[nv ne nd]);
xm = reshape(mean(p,1),[ne nd]);

idx1 = find((xm(:,2) - 2*xm(:,1))>0);
idx2 = find((xm(:,2) - 2*xm(:,1))<0);
t1 = mesh.t(idx1,:);
t2 = mesh.t(idx2,:);

bndexpr = {'all(p(:,2)<min(p0(:,2))+1e-3)', ...
           'all(p(:,2)>max(p0(:,2))-1e-3)',...
           'all(p(:,1)<min(p0(:,1))+1e-3)','true'};     
mesh1 = mkmesh(mesh.p,t1,porder,bndexpr,0,1);

bndexpr = {'all(p(:,2)<min(p0(:,2))+1e-3)', ...
           'all(p(:,2)>max(p0(:,2))-1e-3)',...
           'all(p(:,1)>max(p0(:,1))-1e-3)','true'};     
mesh2 = mkmesh(mesh.p,t2,porder,bndexpr,0,1);

figure(1); clf; meshplot(mesh);
figure(2); clf; meshplot(mesh1);
figure(3); clf; meshplot(mesh2);



