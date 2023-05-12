function mesh = naca0012shock2(poly, poly2, porder)

yu = loginc(linspace(0.04326, 0.7172, 21),0.25);
xu = polyval(poly,yu);

yl = linspace(-0.16794,-0.0597, 5);
xl = polyval(poly2,yl);

du = [0 0.0008 0.002 0.004 0.008 0.016 0.032 0.06];
du = [-du(end:-1:2) du];
theta = -pi/36; nu = [cos(theta) sin(theta)];
[pu,tu] = shockgrid(xu, yu, du, nu, 0);

dl = [0 0.0015 0.0035 0.0075 0.015 0.027];
dl = [-dl(end:-1:2) dl];
theta = pi/150; nl = [cos(theta) sin(theta)];
[pl,tl] = shockgrid(xl, yl, dl, nl, 1);

bndexpr = {'all(p(:,2)<0.05)','all(p(:,2)>0.71)','all(p(:,1)<0.61)','all(p(:,1)>0.61)'};                  
meshl = mkmesh(pu,tu,1,bndexpr,0,1);
p1 = boundarypoints(meshl.p,meshl.f,2);
p2 = boundarypoints(meshl.p,meshl.f,3);
p3 = boundarypoints(meshl.p,meshl.f,4);
pbu = [p3; p1(end-1:-1:1,:); p2(end-1:-1:1,:)];

x1 = loginc(linspace(1.0089304, pbu(1,1), 7),0.6);
y1 = naca12(x1);       
pbu = [x1(:) y1(:); pbu(2:end,:)];
x1 = loginc(linspace(pbu(end,1), 0, 20),1.5);
x1 = [x1(1) 0.530 0.48 x1(4:end)]
x1 = [x1(1:end-1) x1(end-1)/3 x1(end-1)/16 0];
y1 = naca12(x1);       
pbu = [pbu(1:end-1,:); x1(:) y1(:)];

bndexpr = {'all(p(:,2)>-0.065)','all(p(:,2)<-0.16)','all(p(:,1)>0.325)','all(p(:,1)<0.325)'};                  
meshl = mkmesh(pl,tl,1,bndexpr,0,1);
p1 = boundarypoints(meshl.p,meshl.f,2);
p2 = boundarypoints(meshl.p,meshl.f,3);
p3 = boundarypoints(meshl.p,meshl.f,4);
pbl = [p3(end:-1:1,:); p1(2:end,:); p2(2:end,:)];

x1 = loginc(linspace(0, pbl(1,1), 12),1.25);
x1 = [x1(2)/16 x1(2)/3 x1(2:end)]; %x1(end-2:end-1) = [0.24 0.29];
y1 = -naca12(x1);       
pbl = [x1(:) y1(:); pbl(2:end,:)];
x1 = loginc(linspace(pbl(end,1), 1.0089304, 14),0.6);
x1 = x1(1:end-1);
y1 = -naca12(x1);       
pbl = [pbl(1:end-1,:); x1(:) y1(:)];

xd = [-4, 5, -4, 4];
pv1 = [pbu; pbl];
pv2 = [xd(1),xd(3); xd(2),xd(3); xd(2),xd(4); xd(1),xd(4)];
[p, t] = polymesh({pv1,pv2},[1,1],[0,1;1,0],[0.85,1.85]);
idx = find(abs((p(:,1)-0.305).^2 + (p(:,2)+0.06).^2)<1e-6);
p(idx,2) = -naca12(p(idx,1));
% i1 = abs(naca12(p(:,1))-p(:,2)) <1e-4;       
% p(i1,2) = naca12(p(i1,1));
% i1 = abs(naca12(p(:,1))+p(:,2)) <1e-4;       
% p(i1,2) = -naca12(p(i1,1));

% pv = [0.5455 0.711; 0.673 0.6998; 0.3059 -0.1263; 0.3599 -0.1252; 0.3554 -0.1538; 0.9665 6.0e-3; 0.9768 -4.5e-3];
%pv = [0.572 0.162; 0.5426 0.7110; 0.6701 0.6998; 0.6096 0.7167; 0.3281 -0.1680; 0.966 5.9e-3; 0.977 -4.3e-3];
pv = [0.5426 0.7110; 0.6701 0.6998; 0.61 0.717; 0.966 5.9e-3; 0.977 -4.3e-3; 0.442 0.056];
[p,t]=mergeelemementatnode(p, t, pv);

[p,t] = connectmesh(p,t,pu,tu,1e-8);
[p,t] = connectmesh(p,t,pl,tl,1e-8);
[p,t]=fixmesh(p,t);

elemtype = 0;
nodetype = 1;
bndexpr = {'all(sqrt(sum(p.^2,2))<2)','all(sqrt(sum(p.^2,2))>2)'};  

mesh = mkmesh(p,t,porder,bndexpr,elemtype,nodetype);

fb = 1;
mesh.fcurved = (mesh.f(:,end)==-fb);
ic = mesh.fcurved;
mesh.tcurved = false(size(mesh.t,1),1);
mesh.tcurved(mesh.f(ic,end-1)) = true;
mesh.dgnodes = makedgnodes(mesh,@naca12dist,fb);


figure(1); clf; hold on;
simpplot(p,t); axis on;

simpplot(pu,tu); axis on;
simpplot(pl,tl); axis on;
plot(pbl(:,1),pbl(:,2),'o');
plot(pbu(:,1),pbu(:,2),'o');
size(t)
