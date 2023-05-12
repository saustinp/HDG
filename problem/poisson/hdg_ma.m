function [u,q,uhat,v] = hdg_ma(master, mesh, app, u, uhat)
%HDG_POI Solve the poission problem using HDG

[M, Minv, MiB, MiC, MiD, MiE, F, H, G, L, S] = hdg_compute(master, mesh, app);
 
% compute the index mapping
elcon = mesh.elcon;

q = getq(MiB, MiC, u, uhat, elcon);
v = getv(MiD, MiE, u, uhat, elcon);

while (1)
  [AE, FE, Ua, Wa] = assemble(master, mesh, M, MiB, MiC, MiD, MiE, F, H, G, L, S, v, q, u, uhat);

  [duhat, normR] = linearsystem(AE, FE, elcon);
  du = getu(Wa, Ua, duhat, elcon);

  uhat = uhat + duhat;
  u = u + du;
  q = getq(MiB, MiC, u, uhat, elcon);
  v = getv(MiD, MiE, u, uhat, elcon);
  
  normR 
  if normR<1e-8
    break;
  end
end

end

function [AE, FE, Ua, Wa] = assemble(master, mesh, M, MiB, MiC, MiD, MiE, F, H, G, L, S, v, q, u, uhat)

npv = size(M,1);
[npfe, nt] = size(S);
nd = size(MiB,3);
ngv = master.ngv;

shap = squeeze(master.shapvl(:,:,1));
shapmv = master.shapmv(:,:,1);

Qa = zeros(npv*nd,npfe);
Ra = zeros(npv*nd,1);
vg = zeros(ngv, nd, nd);
D = zeros(npv, npv, nd, nd);

Ua = zeros(npv,npfe,nt);
Wa = zeros(npv,nt);
AE = zeros(npfe,npfe,nt);
FE = zeros(npfe,nt);

for i = 1:nt % loop over each element    
    %Mass = M(:,:,i);    
    M_K = H(:,:,i);
    G_K = G(:,:,i);
    I_K = L(:,:,i);
    S_K = S(:,i);
    
    [~, ~, jac] = volgeom(master.shapmv,reshape(mesh.dgnodes(:,:,i),[npv 1 nd]));      
    for l = 1:nd
      vg(:,:,l) = shapmv*v(:,:,l,i);
    end
    % r(w) = (v11*v22 - v12*v21, w)
    Ru = F(:,i) - shap*(master.gwvl.*jac.*(vg(:,1,1).*vg(:,2,2) - vg(:,1,2).*vg(:,2,1)));
    % J = (v22 * dv11, w) + (v11 * dv22, w) - (v21 * dv12, w)  - (v12 * dv21, w)    
    for l = 1:nd
      for m = 1:nd
        D(:,:,l,m)  = shap*diag(master.gwvl.*jac.*vg(:,l,m))*shap';    
      end
    end    
    
    % dV = MiD*dU + MiE*dUhat
    % D11 * dv22 + D22*dv11 - D12*dv21 - D21*dv12
    Di =   D(:,:,1,1)*MiD(:,:,2,2,i) + D(:,:,2,2)*MiD(:,:,1,1,i)  ...
         - D(:,:,1,2)*MiD(:,:,2,1,i) - D(:,:,2,1)*MiD(:,:,1,2,i);
    Ei =   D(:,:,1,1)*MiE(:,:,2,2,i) + D(:,:,2,2)*MiE(:,:,1,1,i)  ...
         - D(:,:,1,2)*MiE(:,:,2,1,i) - D(:,:,2,1)*MiE(:,:,1,2,i);
    Ei = -Ei;
    Ua(:,:,i) = Di\Ei;    
    Wa(:,i) = Di\Ru;    
    
    imap = mesh.elcon(:,i);        
    Rh = S_K - M_K*uhat(imap) - G_K*reshape(q(:,:,i),[npv*nd 1]) + I_K*u(:,i);
    
    for l = 1:nd
      Qa((l-1)*npv+1:l*npv,:) = MiC(:,:,l,i) - MiB(:,:,l,i)*Ua(:,:,i);
      Ra((l-1)*npv+1:l*npv,:) = MiB(:,:,l,i)*Wa(:,i);
    end    
    H_K = M_K + G_K*Qa - I_K*Ua(:,:,i);
    R_K = G_K*Ra + I_K*Wa(:,i);      
    
    AE(:,:,i) = H_K;
    FE(:,i) = R_K + Rh;    
end

end

function [uhat, normR] = linearsystem(AE, FE, elcon)

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

normR = norm(R);

end


