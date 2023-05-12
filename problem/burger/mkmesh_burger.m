S = [0 0; 1 0; 0 1; 1 1];
L = [-1 0; 0 0; -1 1; 0.5 1];
R = [0 0; 1.0 0; 0.5 1; 1 1];
xrect = [-1 1 0 1];

m = 11; n = 11; 
[p,t] = squaremesh(m,n,0,1);
% p(:,1) = xrect(1) + (xrect(2)-xrect(1))*p(:,1);
% p(:,2) = xrect(3) + (xrect(4)-xrect(3))*p(:,2);

p1 = mapp(p,L);
p2 = mapp(p,R);

figure(1);clf; simpplot(p1,t);
figure(2);clf; simpplot(p2,t);

[p,t] = connectmesh(p1,t,p2,t);
bndexpr = {'all(p(:,2)<min(p0(:,2))+1e-3)','all(p(:,1)>max(p0(:,1))-1e-3)', ...
           'all(p(:,2)>max(p0(:,2))-1e-3)','all(p(:,1)<min(p0(:,1))+1e-3)'};     
mesh = mkmesh(p,t,porder,bndexpr,1,1);

