function [u,q,uhat,v,iter] = hdg_ma4(master, mesh, app, u, uhat)
% Extend hdg_ma2 to deal with neumman boundary condition

[M, ~, MiB, MiC, MiD, MiE, F, H, G, L, S, ~, ~, P] = hdg_compute(master, mesh, app);
 
% compute the index mapping
elcon = mesh.elcon;

q = getq(MiB, MiC, u, uhat, elcon);
v = getv(MiD, MiE, u, uhat, elcon);

[AE, ~, Ua, ~, Da] = assemble(M, MiB, MiC, MiD, MiE, F, H, G, L, S);

mesh2 = mkcgmesh(mesh);
mesh2.ib = [];
mesh2.in = 1:size(mesh2.p2,1);

npe = master.npv;
ne = mesh.ne;
[ut, A] = cghelmholtz(master,mesh2,reshape(u,[npe ne]),1e-4); 
il = zeros(npe,npe,ne);
jl = zeros(npe,npe,ne);
for i=1:ne    
    con = mesh2.t2(i,:)';    
    com = repmat(con,[1 npe]);%[con con con con con con];
    il(:,:,i) = com;
    jl(:,:,i) = com';        
end
Fe = zeros(npe,ne);

vold = v;
qold = q;
itermax = 90;
iter = 0;
while (iter<itermax)
  iter = iter + 1;  
  tic
  [FE, Wa] = assembleRHS(master, mesh, app, u, q, v, MiB, Da, G, L, S);
  toc
  [c, uc] = uaverage(Wa, Ua, P, mesh.elcon);  
  tic
  uhat = linearsystem(AE, FE, elcon, app.neumman, uc, c);  
  toc
  u = getu(Wa, Ua, uhat, elcon);
  %sum(u(:).*P(:))
  
  q = getq(MiB, MiC, u, uhat, elcon);
  v = getv(MiD, MiE, u, uhat, elcon);
  
%   for l = 1:master.nd
%     for m = 1:master.nd
%       for i = 1:mesh.ne
%           Fe(:,i) = M(:,:,i)*v(:,l,m,i);          
%       end
%       vl = sparse(reshape(il(:,1,:),npe*ne,1),ones(npe*ne,1),reshape(Fe,npe*ne,1)); 
%       vt = A\vl; vt = full(vt(:));
%       vl = reshape(vt(mesh2.t2'), npe, 1, 1, ne);    
%       v(:,l,m,:) = vl;
% %       figure(3); clf; scaplot(mesh,v(:,1,1,:),[],0,1); axis equal; axis tight; colormap jet;
% %       figure(4); clf; scaplot(mesh,vl(:,1,1,:),[],0,1); axis equal; axis tight; colormap jet;      
% %       pause
%     end    
%   end
  
  if norm(q(:)-qold(:)) <= 1e-6
    break;
  end
  
%   figure(1); clf; scaplot(mesh,u(:,1,:),[],0,1); axis equal; axis tight; colormap jet;
%   figure(2); clf; scaplot(mesh,q(:,1,:),[],0,1); axis equal; axis tight; colormap jet;
%   figure(3); clf; scaplot(mesh,v(:,1,1,:),[],0,1); axis equal; axis tight; colormap jet;
%   figure(4); clf; scaplot(mesh,v(:,2,2,:),[],0,1); axis equal; axis tight; colormap jet;
%   figure(5); clf; scaplot(mesh,v(:,1,2,:),[],0,1); axis equal; axis tight; colormap jet;
%   figure(6); clf; scaplot(mesh,v(:,2,1,:),[],0,1); axis equal; axis tight; colormap jet;
%   detJ = v(:,1,1,:).*v(:,2,2,:) - v(:,1,2,:).*v(:,2,1,:);
%   figure(7); clf; scaplot(mesh,detJ,[-1 1],0,1); axis equal; axis tight; colormap jet;
%   x = mesh.dgnodes(:,1,:); x = x(:);
%   y = mesh.dgnodes(:,2,:); y = y(:);
%   ib1 = find(abs(y+0.5)<1e-6);
%   ib2 = find(abs(y-0.5)<1e-6);
%   ib3 = find(abs(x+0.5)<1e-6);
%   ib4 = find(abs(x-0.5)<1e-6);
%   ib = unique([ib1; ib2; ib3; ib4]);
%   ub = u(:);
%   qx = q(:,1,:); qx = qx(:);
%   qy = q(:,2,:); qy = qy(:);
%   [max(abs((x(ib).^2+y(ib).^2)/2-ub(ib))) max(abs(x(ib)-qx(ib))) max(abs(y(ib)-qy(ib)))]
      
  [iter norm(q(:)-qold(:))]
  vold = v;
  qold = q;  
end

% figure(2); clf; meshplot(mesh,1);    
% mesh2 = mesh; 
% mesh2.dgnodes = q;
% figure(3); clf; meshplot(mesh2,1);


end


function [AE, FE, Ua, Wa, Da] = assemble(M, MiB, MiC, MiD, MiE, F, H, G, L, S)

npv = size(M,1);
[npfe, nt] = size(S);
nd = size(MiB,3);

Qa = zeros(npv*nd,npfe);
Ra = zeros(npv*nd,1);

Da = zeros(npv,npv,nt);
Ua = zeros(npv,npfe,nt);
Wa = zeros(npv,nt);
AE = zeros(npfe,npfe,nt);
FE = zeros(npfe,nt);

for i = 1:nt % loop over each element    
    Mass = M(:,:,i);    
    F_K = F(:,i);
    M_K = H(:,:,i);
    G_K = G(:,:,i);
    I_K = L(:,:,i);
    S_K = S(:,i);
        
    Di = -Mass*MiD(:,:,1,1,i);
    Ei = Mass*MiE(:,:,1,1,i);
    for l = 2:nd
      Di = Di - Mass*MiD(:,:,l,l,i);
      Ei = Ei + Mass*MiE(:,:,l,l,i);
    end        
    % Di*U = Ei*Uhat + F_K
    Da(:,:,i) = Di;
    Ua(:,:,i) = Di\Ei;    
    Wa(:,i) = Di\F_K;    
    
    for l = 1:nd
      Qa((l-1)*npv+1:l*npv,:) = MiC(:,:,l,i) - MiB(:,:,l,i)*Ua(:,:,i);
      Ra((l-1)*npv+1:l*npv,:) = MiB(:,:,l,i)*Wa(:,i);
    end    
    H_K = M_K + G_K*Qa - I_K*Ua(:,:,i);
    R_K = G_K*Ra + I_K*Wa(:,i);      
    
    AE(:,:,i) = H_K;
    FE(:,i) = R_K + S_K;    
end

end

function [FE, Wa] = assembleRHS(master, mesh, app, u, q, v, MiB, Da, G, L, S)

npv = size(Da,1);
[npfe, nt] = size(S);
nd = size(MiB,3);
ngv = master.ngv;

vg = zeros(ngv, nd, nd);
Ra = zeros(npv*nd,1);
Wa = zeros(npv,nt);
FE = zeros(npfe,nt);

source   = str2func(app.source);

shap = squeeze(master.shapvl(:,:,1));
for i = 1:nt % loop over each element    
    G_K = G(:,:,i);
    I_K = L(:,:,i);
    S_K = S(:,i);
            
    % form F_K
    [xg, ~, jac] = volgeom(master.shapmv,reshape(mesh.dgnodes(:,:,i),[npv 1 nd]));      
    ug = master.shapmv(:,:,1)*u(:,1,i);
    qg = master.shapmv(:,:,1)*q(:,:,i);
    sg = source(xg, [ug qg], app.param);        
    for l = 1:nd
      vg(:,:,l) = master.shapmv(:,:,1)*v(:,:,l,i);
    end        
    fg = -sqrt(2*sg + vg(:,1,1).^2 + vg(:,2,2).^2 + vg(:,1,2).^2 + vg(:,2,1).^2);
    F_K = shap*(master.gwvl.*jac.*fg);    
    
    Wa(:,i) = Da(:,:,i)\F_K;        
    for l = 1:nd
      Ra((l-1)*npv+1:l*npv,:) = MiB(:,:,l,i)*Wa(:,i);
    end    
    R_K = G_K*Ra + I_K*Wa(:,i);      
    
    FE(:,i) = R_K + S_K;    
end

end


function uhat = linearsystem(AE, FE, elcon, neumann, uc, c)

[npfe, nt] = size(FE);

il = zeros(npfe,npfe,nt);
jl = zeros(npfe,npfe,nt);
for i=1:nt
    con = repmat((elcon(:,i)'-1),1,1)+ones(1,npfe);
    con = reshape(con,npfe,1);
    il(:,:,i) = repmat(con ,1,npfe);
    jl(:,:,i) = repmat(con',npfe,1);        
end
H = sparse(reshape(il,(npfe)^2*nt,1),reshape(jl,(npfe)^2*nt,1),reshape(AE,(npfe)^2*nt,1));
R = sparse(reshape(il(:,1,:),(npfe)*nt,1),ones((npfe)*nt,1),reshape(FE,(npfe)*nt,1));                    

if neumann==1
%   n = length(R);
%   e = zeros(1,n); 
%   e(1) = 1;
%   e = sparse(e);  
%   R(end+1) = 0;
%   H(end+1,:) = e;    
%   uhat = H\R;  
  
  R(end+1) = -c;
  H(end+1,:) = uc;  
  H(:,end+1) = [uc(:); 0*max(abs(uc(:)))];
  ut = H\R;
  uhat = ut(1:end-1);    
  %[ut(1) ut(end)]
  
%   H(:,1) = [];
%   ut = H\R;
%   uhat = [0; ut];

%     H(1,2:end) = 0;
%     %H(2:end,1) = 0;
%     R(1) = 0;   
%     uhat = H\R;    
else  
  uhat = H\R;
end

end


