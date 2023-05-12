function [AE, FE, We] = hdg_assembly(master, mesh, app, AE, MiB, MiC, De, Ue, F, H, G, L, S)

npv = size(Da,1);
[npfe, nt] = size(S);
nd = size(MiB,3);
nfe = size(master.perm,2);

Qa = zeros(npv*nd,npfe);
Ra = zeros(npv*nd,1);
We = zeros(npv,nt);
FE = zeros(npfe,nt);

for i = 1:nt % loop over each element    
    F_K = F(:,i);
    G_K = G(:,:,i);
    I_K = L(:,:,i);
    S_K = S(:,i);    
                    
    We(:,i) = De(:,:,i)\F_K;        
    for l = 1:nd
      Ra((l-1)*npv+1:l*npv,:) = MiB(:,:,l,i)*We(:,i);
    end    
    
    nbc = 0; % check boundary conditions
    for j = 1:nfe 
      ibf = -mesh.f(abs(mesh.t2f(i,j)),end);
      if (ibf>0) && (app.bcm(ibf) > 1)          
        nbc = 1;        
      end
    end
    
    if nbc == 1      
      for l=1:nd      
        Qa((l-1)*npv+1:l*npv,:) = MiC(:,:,l,i) - MiB(:,:,l,i)*Ue(:,:,i);
      end
      AE(:,:,i) = H(:,:,i) + G_K*Qa - I_K*Ue(:,:,i);           
    end        
    
    FE(:,i) = G_K*Ra + I_K*We(:,i) + S_K;    
end

end
