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
tau = 1;
av = 0.0;

app.arg = {kappa,c,av,1/(ngrid*porder),tau};
app.bcm = [4;1;2;1];
app.bcs = [0;0;0;1]; %[1,0,0;1,0,0];
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
mkmesh_burger2;
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);

% initial solution
UDG = initu(mesh,{0;0;0});
UH=inituhat(master,mesh.elcon,UDG,1);

% CG discretization for AV field
[cgnodes,cgelcon,rowent2elem,colent2elem,cgent2dgent] = mkcgent2dgent(mesh.dgnodes(:,1:2,:),1e-8);
R = 0.1; wR = 0.1;
[nodeR,weightR] = nodelistcg(mesh,cgnodes,cgelcon,R,wR);

% HDG solver for constant viscosity field
mesh.dgnodes(:,3,:) = 0.1;
[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);
figure(1); clf;scaplot(mesh,UDG(:,1,:),[],2,1); axis equal; axis tight;

% % HDG solver for constant viscosity field
mesh.dgnodes(:,3,:) = 0.01;
[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);
figure(1); clf;scaplot(mesh,UDG(:,1,:),[],2,1); axis equal; axis tight;

mesh.dgnodes(:,3,:) = 0.002;
[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);
figure(1); clf;scaplot(mesh,UDG(:,1,:),[],2,1); axis equal; axis tight;

% mesh.dgnodes(:,3,:) = 0.001;
% [UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);
% figure(1); clf;scaplot(mesh,UDG(:,1,:),[],2,1); axis equal; axis tight;

% mesh.dgnodes(:,3,:) = 0;
% [UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);
% figure(1); clf;scaplot(mesh,UDG(:,1,:),[],2,1); axis equal; axis tight;


% 
% mesh.dgnodes(:,3,:) = 0.00001;
% [UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);
% figure(1); clf;scaplot(mesh,UDG(:,1,:),[],2,1); axis equal; axis tight;

% alpha = 100;
% av = [0.1 0.09 0.08];
% for i = 1:length(av)
%     div = app.arg{4}*(UDG(:,2,:) + UDG(:,3,:));
%     a = avfield(div, cgelcon, cgent2dgent, rowent2elem, nodeR, weightR, alpha);
%     mesh.dgnodes(:,3,:) = av(i)*a;
%     figure(2); clf;scaplot(mesh, mesh.dgnodes(:,3,:) , [],2,1); axis equal; axis tight;
%         
%     s = ucg2udg(udg2ucg(div, cgent2dgent, rowent2elem), cgelcon);
%     figure(3); clf;scaplot(mesh, s , [0 0.5],2,1); axis equal; axis tight;
%     
%     [UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);
%     figure(4); clf;scaplot(mesh,UDG(:,1,:),[],2,1); axis equal; axis tight;        
% end
% 
% [x,y] = meshgrid(-0.2:0.0025:0.8, 0:0.05:1);
% [pn, un] = shocklocation(mesh, s, x, y, 4);
% poly = polyfit(pn(:,1) ,pn(:,2), 6);
% x1 = linspace(pn(1,1),pn(end,1),1000);
% y1 = polyval(poly,x1);
% 
% figure(1); clf;scaplot(mesh,UDG(:,1,:),[],2,1); axis equal; axis tight;
% hold on;
% plot(pn(:,1), pn(:,2), 'or', 'LineWidth', 2);
% plot(x1, y1, '-r', 'LineWidth', 2); 
% 
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
% 
