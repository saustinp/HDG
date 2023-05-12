function [UDG, UH, ACG, mine, minf, ming, lambda] = avstart(master, mesh, app, UDG0, UH0, avfield, S0, kappa0)

alpha = 100;
href = 1;
nd = mesh.nd;
porder = mesh.porder;

if nd == 1
    np1 = porder;
    A = tensorproduct(master.plocvl, mesh.porder);    
elseif nd==2
    if mesh.elemtype==0
        np1 = porder*(porder+1)/2;
        A = koornwinder(master.plocvl, mesh.porder);    
    else
        np1 = porder*porder;
        A = tensorproduct(master.plocvl, mesh.porder);
    end
else
    if mesh.elemtype==0
        np1 = porder*(porder+1)*(porder+2)/6;
        A = koornwinder(master.plocvl, mesh.porder);    
    else
        np1 = porder*porder*porder;
        A = tensorproduct(master.plocvl, mesh.porder);
    end
end
he = meshsize(mesh);

div = divergence(UDG0, href);
sm = squeeze(mean(div,1));
idx = find(sm > S0); % shock elements
hs = he(:,idx);
hk = (kappa0*min(hs(:)))^2;

% solve the Helmholtz equation
%s = (div).*(atan(alpha*(div))/pi + 0.5) - atan(alpha)/pi + 0.5;    
divmax = max(div(:));
s = limiting(div,0,divmax/2,alpha,0);
s = cgpoisson(mesh, master, s, [hk 1.0]);        
a = (s-S0).*(atan(alpha*(s-S0))/pi + 0.5) - atan(alpha)/pi + 0.5;    
a = reshape(a(mesh.t2'), master.npe, 1, mesh.ne);
a = a/max(a(:));

% s = (div-S0).*(atan(alpha*(div-S0))/pi + 0.5) - atan(alpha)/pi + 0.5;    
% a = cgpoisson(mesh, master, s, [hk 1.0]);        
% a = reshape(a(mesh.t2'), master.npe, 1, mesh.ne);
% figure(1); clf; scaplot(mesh, a,[],2,0); axis on;
% figure(2); clf; scaplot(mesh, a2,[],2,0); axis on;
% figure(3); clf; scaplot(mesh, a-a2,[],2,0); axis on;
% pause
% compute lambda
% c1 = calnorm(mesh, master, avfield, idx);
% [c2, sumK] = calnorm(mesh, master, a, idx);
% lambda = c1/c2;
% c3 = calnorm(mesh, master, lambda*a, idx);
% [c1 c2 c3 lambda c1/sumK]
c1 = avfield(:,1,idx);
c2 = a(:,1,idx);
lambda = max(c1(:))/max(c2(:));

% solve the conservation laws
ACG = lambda*a; 
mesh.dgnodes(:,nd+1,:) = ACG;

figure(1); clf; scaplot(mesh, avfield,[],2,0); axis on;
figure(2); clf; scaplot(mesh, ACG,[],2,0); axis on;
figure(3); clf; scaplot(mesh, a,[],2,0); axis on;
figure(4); clf; scaplot(mesh, div,[],2,0); axis on;
[UDG, UH] = hdg_solve(master,mesh,app,UDG0,UH0,[]);        

pres = UDG(:,1,:);
u = A\squeeze(pres);
u((np1+1):end,:) = 0;
U = reshape(A*u,[master.npe 1 mesh.ne]);    
err = calerror(mesh,master,pres,U,1);
mine = max(err(:));

pres = eulereval(UDG,'p',app.arg{1},app.arg{2});
u = A\squeeze(pres);
u((np1+1):end,:) = 0;
U = reshape(A*u,[master.npe 1 mesh.ne]);    
err = calerror(mesh,master,pres,U,1);
minf = max(err(:));    

pres = eulereval(UDG,'M',app.arg{1},app.arg{2});
u = A\squeeze(pres);
u((np1+1):end,:) = 0;
U = reshape(A*u,[master.npe 1 mesh.ne]);    
err = calerror(mesh,master,pres,U,1);
ming = max(err(:));    

[mine minf ming]       
    
function [err, sumK] = calnorm(mesh, master, UDG, idx)

npv = size(UDG,1);
ngv = master.ngv;
nd  = master.nd;
shapvt    = squeeze(master.shapvt(:,:,1));
dshapvt   = reshape(permute(master.shapvt(:,:,2:nd+1),[1 3 2]),[ngv*nd npv]);

err = 0;
sumK = 0;
for j = 1:length(idx)
    i = idx(j);
    dg = mesh.dgnodes(:,:,i);
    
    % compute the Jacobian matrix at Gauss points: dx/dxi
    Jg = dshapvt*dg(:,1:nd);
    Jg = reshape(Jg,[ngv nd nd]);        
    jac = volgeom(Jg);                   
    udgg = shapvt*UDG(:,1,i);    
    err = err + (master.gwvl.*jac)'*udgg;      
    sumK = sumK + sum(master.gwvl.*jac);
end


function err = calerror(mesh,master,UDG, U0, p)

[npv, nc, ne] = size(UDG);
ngv = master.ngv;
nd  = master.nd;

shapvt    = squeeze(master.shapvt(:,:,1));
dshapvt   = reshape(permute(master.shapvt(:,:,2:nd+1),[1 3 2]),[ngv*nd npv]);

err = zeros(nc,ne);
for i = 1:ne
    dg = mesh.dgnodes(:,:,i);
    
    % compute the Jacobian matrix at Gauss points: dx/dxi
    Jg = dshapvt*dg(:,1:nd);
    Jg = reshape(Jg,[ngv nd nd]);        
    jac = volgeom(Jg);                   
    udgg = shapvt*UDG(:,:,i);    
    udge = shapvt*U0(:,:,i);        
    for j = 1:nc
        err(j,i) = (master.gwvl.*jac)'*(abs(udgg(:,j)-udge(:,j)).^p)/sum(master.gwvl.*jac);  
    end    
end
%err  = sqrt(err);

function [jac] = volgeom(Jg)

nd  = size(Jg,2);
switch nd
    case 1
        jac = Jg;        
    case 2
        jac = Jg(:,1,1).*Jg(:,2,2) - Jg(:,1,2).*Jg(:,2,1);        
    case 3
        jac = Jg(:,1,1).*Jg(:,2,2).*Jg(:,3,3) - Jg(:,1,1).*Jg(:,3,2).*Jg(:,2,3)+ ...
              Jg(:,2,1).*Jg(:,3,2).*Jg(:,1,3) - Jg(:,2,1).*Jg(:,1,2).*Jg(:,3,3)+ ...
              Jg(:,3,1).*Jg(:,1,2).*Jg(:,2,3) - Jg(:,3,1).*Jg(:,2,2).*Jg(:,1,3);                    
    otherwise
        error('Dimension is not implemented');
end





