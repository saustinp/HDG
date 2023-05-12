setapplicationpath('FM/burgersav');

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
mesh = mkmesh_unstructuredrect(0.25,porder);
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);

% initial solution
UDG = initu(mesh,{0;0;0});
UH=inituhat(master,mesh.elcon,UDG,1);

% CG discretization for AV field
[cgnodes,cgelcon,rowent2elem,colent2elem,cgent2dgent] = mkcgent2dgent(mesh.dgnodes(:,1:2,:),1e-8);
R = 0.25; wR = 0.1;
[nodeR,weightR] = nodelistcg(mesh,cgnodes,cgelcon,R,wR);

mesh.dgnodes(:,3,:) = 0.01;
[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);
figure(1); clf;scaplot(mesh,UDG(:,1,:),[],2,1); axis equal; axis tight;

mesh.dgnodes(:,3,:) = 0.0;
alpha = 100;
av = [0.2];
for i = 1:length(av)
    div = app.arg{4}*(UDG(:,2,:) + UDG(:,3,:));
    a = avfield(div, cgelcon, cgent2dgent, rowent2elem, nodeR, weightR, alpha);
    mesh.dgnodes(:,3,:) = av(i)*a;
    figure(2); clf;scaplot(mesh, mesh.dgnodes(:,3,:) , [],2,1); axis equal; axis tight;
        
    s = ucg2udg(udg2ucg(div, cgent2dgent, rowent2elem), cgelcon);
    figure(3); clf;scaplot(mesh, s , [0 0.5],2,1); axis equal; axis tight;
    
    [UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);
    figure(4); clf;scaplot(mesh,UDG(:,1,:),[],2,1); axis equal; axis tight;        
end

% a = [1 1/2 1/4]*0.01;
% for i = 1:length(a)
%     mesh.dgnodes(:,3,:) = a(i);
%     [UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);
%     figure(1); clf;scaplot(mesh,UDG(:,1,:),[],2,2); axis equal; axis tight;
% end
% 

% mesh.dgnodes(:,3,:) = 0.01/16;
% [UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);
% figure(1); clf;scaplot(mesh,UDG(:,1,:),[],2,1); axis equal; axis tight;
% 
% mesh.dgnodes(:,3,:) = 0.01/4;
% [UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);
% figure(1); clf;scaplot(mesh,UDG(:,1,:),[],2,1); axis equal; axis tight;
% 
% mesh.dgnodes(:,3,:) = 0.01;
% [UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);
% figure(1); clf;scaplot(mesh,UDG(:,1,:),[],2,1); axis equal; axis tight;
% 
% div = app.arg{4}*(UDG(:,2,:) + UDG(:,3,:));
% a = avfield(div, cgelcon, cgent2dgent, rowent2elem, nodeR, weightR, alpha);
% mesh.dgnodes(:,3,:) = 0.1*a;
% figure(2); clf;scaplot(mesh, mesh.dgnodes(:,3,:) , [],2,1); axis equal; axis tight;
% 
% s = ucg2udg(udg2ucg(div, cgent2dgent, rowent2elem), cgelcon);
% figure(3); clf;scaplot(mesh, s , [0 0.5],2,1); axis equal; axis tight;
% 
% [UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);
% 
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
