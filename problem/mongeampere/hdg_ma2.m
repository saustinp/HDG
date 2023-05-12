function [u,q,uhat,v,iter] = hdg_ma2(master, mesh, app, u, uhat)
%HDG_POI Solve the poission problem using HDG

[M, Minv, MiB, MiC, MiD, MiE, F, H, G, L, S] = hdg_compute(master, mesh, app);
 
% compute the index mapping
elcon = mesh.elcon;

q = getq(MiB, MiC, u, uhat, elcon);
v = getv(MiD, MiE, u, uhat, elcon);

[AE, ~, Ua, ~, Da] = assemble(M, MiB, MiC, MiD, MiE, F, H, G, L, S);

vold = v;
qold = q;
itermax = 90;
iter = 0;
while (iter<itermax)
  iter = iter + 1;
  
  [FE, Wa] = assembleRHS(master, mesh, app, v, MiB, Da, G, L, S);

  uhat = linearsystem(AE, FE, elcon);

  u = getu(Wa, Ua, uhat, elcon);
  q = getq(MiB, MiC, u, uhat, elcon);
  v = getv(MiD, MiE, u, uhat, elcon);
  if norm(q(:)-qold(:)) <= 1e-8
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
%   pause
  
  [iter norm(q(:)-qold(:))]
  vold = v;  
  qold = q;
end

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

function [FE, Wa] = assembleRHS(master, mesh, app, v, MiB, Da, G, L, S)

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
    sg = source(xg, [xg xg(:,1)], app.param);    
    for l = 1:nd
      vg(:,:,l) = master.shapmv(:,:,1)*v(:,:,l,i);
    end        
    fg = -sqrt(-2*sg + vg(:,1,1).^2 + vg(:,2,2).^2 + vg(:,1,2).^2 + vg(:,2,1).^2);
    F_K = shap*(master.gwvl.*jac.*fg);    
    
    Wa(:,i) = Da(:,:,i)\F_K;        
    for l = 1:nd
      Ra((l-1)*npv+1:l*npv,:) = MiB(:,:,l,i)*Wa(:,i);
    end    
    R_K = G_K*Ra + I_K*Wa(:,i);      
    
    FE(:,i) = R_K + S_K;    
end

end


function uhat = linearsystem(AE, FE, elcon)

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
uhat = H\R;

end


