%setapplicationpath('FM/burgersav');

porder = 4;
nstage = 1;
torder = 1;
ngrid = 12;
hybrid = 'hdg';

clear app;
app.source = 'source';
app.flux = 'flux';
app.fbou = 'fbou';
app.fhat = 'fhat';
app.hybrid = hybrid;
app.localsolve=1;
ntime  = 5;
dt = linspace(0.01,10,ntime);

kappa = 0.0;
c = [0.5,0.0];
tau = 0.5;
av = 0.0;

app.arg = {kappa,c,av,1/(ngrid*porder),tau};
app.bcm = [4;1;2;1];
app.bcs = [0;0;0;0]; %[1,0,0;1,0,0];
app.bcd  = [1,1,1,1];  % 2: Slip wall, 1: Far-field
app.bcv  = [0; 0; 0; 0];

app.tdep = false;
app.wave = false;
app.alag = false;
app.flg_q = 1;
app.flg_p = 0;
app.flg_g = 0;

app.nd = 2;
app.nch  = 1;                       % Number of componets of UH
app.nc   = app.nch*(app.nd+1);    % Number of componeents of UDG
app.ncu  = 1;

app.time = [];
app.dtfc = [];
app.alpha = [];
app.itmax = 1;
app.fc_q = 1;
app.fc_p = 0;
app.fc_u = 0;
app.tdep = false;
app.nc = 3;

% mesh and master
%mesh = mkmesh_rect(ngrid*2+1,ngrid+1,porder,0,[-1 1 0 1],1,1);
mesh = mkmesh_shockrect(porder);
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);

% initial solution
UDG = initu(mesh,{0.0;0;0});
UH=inituhat(master,mesh.elcon,UDG,1);

% CG discretization for AV field
[cgnodes,cgelcon,rowent2elem,colent2elem,cgent2dgent] = mkcgent2dgent(mesh.dgnodes(:,1:2,:),1e-8);
R = 0.1; wR = 0.1;
[nodeR,weightR] = nodelistcg(mesh,cgnodes,cgelcon,R,wR);

% HDG solver for constant viscosity field
mesh.dgnodes(:,3,:) = 0.01;
[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);
figure(1); clf;scaplot(mesh,UDG(:,1,:),[],2,0); axis equal; axis tight;

A = tensorproduct(master.plocvl,porder);

beta = 0;
alpha = 100;

mesh=mkcgmesh(mesh);
div = UDG(:,2,:) + beta*UDG(:,3,:);
s = ucg2udg(udg2ucg(div, cgent2dgent, rowent2elem), cgelcon);
mesh.ib = [];
mesh.in = 1:size(mesh.p2,1);
u = cgpoisson(mesh,master,s,[0.0001 1.0]);
[ne, npe] = size(mesh.t2);
u = u(mesh.t2');
u = reshape(u,npe,1,ne);
figure(2); clf;scaplot(mesh,u,[0 20],2,1); axis equal; axis tight;
pause

av = [0.05 0.025 0.0125 0.0125/2]/2;
for i = 1:length(av)
    div = app.arg{4}*(UDG(:,2,:) + beta*UDG(:,3,:));
    a = avfield(div, cgelcon, cgent2dgent, rowent2elem, nodeR, weightR, alpha);
    mesh.dgnodes(:,3,:) = av(i)*a;
    figure(1); clf;scaplot(mesh, mesh.dgnodes(:,3,:) , [],2,0); axis equal; axis tight;
        
    s = ucg2udg(udg2ucg(div, cgent2dgent, rowent2elem), cgelcon);
    figure(2); clf;scaplot(mesh, s , [0 0.5],2,0); axis equal; axis tight;
    
    [UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);
    UDGprev = UDG;
    UHprev = UH;
    figure(3); clf;scaplot(mesh,UDG(:,1,:),[],2,0); axis equal; axis tight;       
    
%     [~,~,e1,e2] = dgprojection2(master,mesh,UDG(:,1,:),porder-1);
%     [max(e1) max(e2)]
    
    u=A\squeeze(UDG(:,1,:));
    u1=u(1:(porder)*(porder),:);
    e1=sum(abs(u1),1);
    e=sum(abs(u),1);
    
%     u2 = 0*u;
%     u2(1:(porder)*(porder),:) = u1;
%     UDG1 = reshape(A*u2,[master.npe 1 mesh.ne]);
%    
%     figure(4); clf;scaplot(mesh,UDG1(:,1,:),[],2,0); axis equal; axis tight;       
%     [UDG(:,1,1) UDG1(:,1,1)]
    
    [i av(i) min(e1./e)]
%     if min(e1./e) < 0.5
%         break;
%     end
end
UDG = UDGprev;
UH = UHprev;
pause

div = UDG(:,2,:) + beta*UDG(:,3,:);
s = ucg2udg(udg2ucg(div, cgent2dgent, rowent2elem), cgelcon);
figure(3); clf;scaplot(mesh, s , [-1 10],2,1); axis equal; axis tight;

s0 = 2; s1 = 4; nref = 2; d = 2; n = 20; m = 1000; polydegree=6;
[px, sx, idx] = shockpoints(mesh, s, s0, s1, nref);
[pp, pl, pr, pa] = shockcells(px, d, n, m);
x = squeeze(pa(:,1,:))';
y = squeeze(pa(:,2,:))';
[pn, un] = shocklocation(mesh, s, idx, x, y, 4);
poly = polyfit(pn(:,1), pn(:,2), polydegree);

k = 12; v = [ones(k,1)  zeros(k,1)];
x1 = linspace(pn(1,1),pn(end,1),k);
y1 = polyval(poly,x1);
y1 = max(y1,0);
y1 = min(y1,1);
ds = [-0.04 -0.02 -0.01 -0.005 0 0.005 0.01 0.02 0.04];
sx = zeros(k, length(ds));
sy = 0*sx;
for i = 1:length(ds)
    sx(:,i) = x1(:) + ds(i)*v(:,1);
    sy(:,i) = y1(:) + ds(i)*v(:,2);
end

figure(1); clf;scaplot(mesh,UDG(:,1,:),[],2,1); axis equal; axis tight;
hold on;
plot(pn(:,1), pn(:,2), 'or', 'LineWidth', 2);
plot(x1, y1, '-r', 'LineWidth', 2); 
plot(sx, sy, '*g');


% div = app.arg{4}*(UDG(:,2,:) + beta*UDG(:,3,:));
% s = avfield(div, cgelcon, cgent2dgent, rowent2elem, nodeR, weightR, alpha);
% figure(1); clf;scaplot(mesh, s , [0 0.08],2,1); axis equal; axis tight;
% 
% s0 = 0.02; s1 = 0.05; nref = 3;
% n = 40; k = [100 3]; 
% [px, sx, idx] = shockpoints(mesh, s, s0, s1, nref);
% [pp, pl, pr, pa] = boundingboxes(px, 2, n, k);
% 
% x = reshape(permute(pa,[1 2 4 3]), [size(pa,1)*size(pa,2)*size(pa,4) size(pa,3)]);
% [udgx] = dgfieldatx(mesh, s, idx, x, 6);
% udgx = reshape(udgx,[k(1) k(2) size(pa,4)]);
% [a,b] = max(udgx,[],1);
% a = squeeze(a);
% b = squeeze(b);
% c = zeros(k(2), n, 2);
% for i = 1:n
%     for j = 1:k(2)
%         c(j,i,:) = pa(b(j,i),j,:,i);
%     end
% end
% figure(2);clf;plot(c(:,:,1),c(:,:,2),'o')
% 
% figure(3); clf;scaplot(mesh, s , [],2,1); axis equal; axis tight;
% hold on;
% plot(px(:,1), px(:,2), 'or', 'LineWidth', 2);
% for i = 1:size(pp,1)
%     plot([pp(i,1) pp(i,2) pp(i,2) pp(i,1) pp(i,1)],[pp(i,3) pp(i,3) pp(i,4) pp(i,4) pp(i,3)],'-b','LineWidth',1);
% end
% for i = 1:size(pa,4)
%     plot(pa(:,:,1,i),pa(:,:,2,i),'og','LineWidth',3);
% end
% 
% figure(2);clf;
% hold on;
% for i = 2:size(pa,4)
%     p1 = pa(:,:,:,i);
%     u1 = udgx(:,:,i);
%     for j = 1:k(2)
%         plot3(p1(:,j,1),p1(:,j,2),u1(:,j));
%     end
% end
% 
% 
% div = UDG(:,2,:) + beta*UDG(:,3,:);
% s = ucg2udg(udg2ucg(div, cgent2dgent, rowent2elem), cgelcon);
% figure(3); clf;scaplot(mesh, s , [-1 10],2,1); axis equal; axis tight;
% 
% [px, sx, idx] = shockpoints(mesh, s, 2, 4, 2);
% 
% [x,y] = meshgrid(-0.2:0.0025:0.8, 0:0.05:1);
% [pp, pl, pr, pa] = shockcells(px, 2, 20, 200);
% x = squeeze(pa(:,1,:))';
% y = squeeze(pa(:,2,:))';
% [pn, un] = shocklocation(mesh, s, x, y, 4);
% poly = polyfit(pn(:,1) ,pn(:,2), 6);
% x1 = linspace(pn(1,1),pn(end,1),1000);
% y1 = polyval(poly,x1);
% 
% figure(1); clf;scaplot(mesh,UDG(:,1,:),[],2,1); axis equal; axis tight;
% hold on;
% plot(pn(:,1), pn(:,2), 'or', 'LineWidth', 2);
% plot(x1, y1, '-r', 'LineWidth', 2); 

% x1 = linspace(pn(1,1),pn(end,1),10);
% y1 = polyval(poly,x1);
% xy1 = [x1(:) y1(:)];
% pv1 = [-1 0; xy1; -1 1];
% [p1,t1]=polymesh({pv1},[1],[0,1],[0.2,1.3]);
% figure(2); clf; simpplot(p1,t1);
% 
% pv2 = [1 0; 1 1; xy1(end:-1:1,:)];
% [p2,t2]=polymesh({pv2},[1],[0,1],[0.2,1.3]);
% figure(3); clf; simpplot(p2,t2);
% hold on; simpplot(p1,t1);


% [pn, un] = shockcurve(mesh, s, 0.15, 4);
% figure(3); clf;scaplot(mesh, s , [0 0.15],2,1); axis equal; axis tight;
% hold on;
% plot(pn(:,1), pn(:,2), 'or', 'LineWidth', 2);
% 
% figure(2); clf;scaplot(mesh, mesh.dgnodes(:,3,:) , [],2,1); axis equal; axis tight;
% hold on;
% plot(pn(:,1), pn(:,2), 'or', 'LineWidth', 2);
% 
% figure(1); clf;scaplot(mesh,UDG(:,1,:),[],2,1); axis equal; axis tight;
% hold on;
% plot(pn(:,1), pn(:,2), 'or', 'LineWidth', 2);
% 
% grad = sqrt(UDG(:,2,:).^2 + UDG(:,3,:).^2);
% s = ucg2udg(udg2ucg(grad, cgent2dgent, rowent2elem), cgelcon);
% [pn, un] = shockcurve(mesh, s, 15, 4);
% figure(3); clf;scaplot(mesh, s, [0 30],2,1); axis equal; axis tight;
% hold on; plot(pn(:,1), pn(:,2), 'or', 'LineWidth', 2);
% 
% 

% [UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);
% figure(2); clf;scaplot(mesh,UDG(:,1,:),[],2,1); axis equal; axis tight;

% av = logdec(linspace(1e-3,1e-2,10),3);
% for i = 1:length(av)       
%     kappa = kappa/(3^i);
%     app.arg = {kappa,c,av(i),1/(ngrid*porder),tau};
%     [UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);
%     figure(1); clf;scaplot(mesh,UDG(:,1,:),[],2,1); axis equal; axis tight;    
% 
%     s = app.arg{4}*(UDG(:,2,:) + beta*UDG(:,3,:));
%     a = s.*(atan(alpha*s)/pi + 0.5) - atan(alpha)/pi + 0.5;
%     figure(2); clf;scaplot(mesh,a,[0 0.4],2,1); axis equal; axis tight;
%     
%     s = dg2cg2(app.arg{4}*(UDG(:,2,:) + beta*UDG(:,3,:)), cgelcon, colent2elem, rowent2elem);
%     a = s.*(atan(alpha*s)/pi + 0.5) - atan(alpha)/pi + 0.5;
%     figure(3); clf;scaplot(mesh,a,[],2,1); axis equal; axis tight;
% 
%     s = dg2cg(app.arg{4}*(UDG(:,2,:) + beta*UDG(:,3,:)), cgelcon, cgent2dgent, rowent2elem);
%     a = s.*(atan(alpha*s)/pi + 0.5) - atan(alpha)/pi + 0.5;
%     figure(4); clf;scaplot(mesh,a,[0 0.4],2,1); axis equal; axis tight;
%     
%     [i av(i)]
% end
% 
% s = dg2cg(app.arg{4}*(UDG(:,2,:) + beta*UDG(:,3,:)), cgelcon, cgent2dgent, rowent2elem);
% S = s;
% a = s.*(atan(alpha*s)/pi + 0.5) - atan(alpha)/pi + 0.5;
% figure(4); clf;scaplot(mesh,a,[0 0.4],2,1); axis equal; axis tight;
% 
% scg = udg2ucg(app.arg{4}*(UDG(:,2,:) + beta*UDG(:,3,:)), cgent2dgent, rowent2elem);
% s = filtering(scg, nodeR, weightR);
% s = ucg2udg(s, cgelcon);
% a = s.*(atan(alpha*s)/pi + 0.5) - atan(alpha)/pi + 0.5;
% figure(5); clf;scaplot(mesh,a,[],2,1); axis equal; axis tight;
% 
% 
% elemlist = elementlist(mesh, R*2);
% s = filtering(sdg, mesh.dgnodes, elemlist, R, wR);
% a = s.*(atan(alpha*s)/pi + 0.5) - atan(alpha)/pi + 0.5;
% figure(5); clf;scaplot(mesh,a,[],2,1); axis equal; axis tight;
% 
% korder = 1;
% mesh1 = mkmesh_rect(ngrid*2+1,ngrid+1,korder,0,[-1 1 0 1],1,1);
% master1 = mkmaster(mesh1,2*korder);
% [master1,mesh1] = preprocess(master1,mesh1,hybrid);
% s = dgprojection(master1,mesh1,S,porder);
% a = s.*(atan(alpha*s)/pi + 0.5) - atan(alpha)/pi + 0.5;
% figure(5); clf;scaplot(mesh1,a,[0 0.4],2,1); axis equal; axis tight;

% s = UDG(:,2,:) + beta*UDG(:,3,:);
% alpha = 100;
% a = s.*(atan(alpha*s)/pi + 0.5) - atan(alpha)/pi + 0.5;
% figure(3); clf;scaplot(mesh,a,[],2,1); axis equal; axis tight;
% 
% [cgnodes,cgelcon,rowent2elem,colent2elem,cgent2dgent] = mkcgent2dgent(mesh.dgnodes,1e-8);
% 
% UCG = dg2cg2(UDG, cgelcon, colent2elem, rowent2elem);
% s = UCG(:,2,:) + UCG(:,3,:);
% a = s.*(atan(alpha*s)/pi + 0.5) - atan(alpha)/pi + 0.5;
% figure(4); clf;scaplot(mesh,a,[],2,1); axis equal; axis tight;
% 
% [UCG,vcg] = dg2cg(UDG, cgelcon, cgent2dgent, rowent2elem);
% s = UCG(:,2,:) + UCG(:,3,:);
% a = s.*(atan(alpha*s)/pi + 0.5) - atan(alpha)/pi + 0.5;
% figure(5); clf;scaplot(mesh,a,[],2,1); axis equal; axis tight;
% 
% s = dg2cg2(UDG(:,2,:) + beta*UDG(:,3,:), cgelcon, colent2elem, rowent2elem);
% a = s.*(atan(alpha*s)/pi + 0.5) - atan(alpha)/pi + 0.5;
% figure(4); clf;scaplot(mesh,a,[],2,1); axis equal; axis tight;
% 
% s = dg2cg(UDG(:,2,:) + beta*UDG(:,3,:), cgelcon, cgent2dgent, rowent2elem);
% a = s.*(atan(alpha*s)/pi + 0.5) - atan(alpha)/pi + 0.5;
% figure(5); clf;scaplot(mesh,a,[],2,1); axis equal; axis tight;


% % HDG postprocessing 
% mesh1 = mkmesh_square(ngrid,ngrid,porder+1);
% master1 = mkmaster(mesh1,2*(porder+1));
% [master1,mesh1] = preprocess(master1,mesh1,hybrid);
% UDGstar = postprocessnd(master,mesh,master1,mesh1,UDG);
% 
% figure(1); clf;scaplot(mesh,UDG(:,1,:),[],2,1); axis equal; axis tight;
% figure(2); clf;scaplot(mesh1,UDGstar(:,1,:),[],2,1); axis equal; axis tight;
