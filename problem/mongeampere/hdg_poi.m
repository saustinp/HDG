function [u,q,uhat,v] = hdg_poi(master, mesh, app)
%HDG_POI Solve the poission problem using HDG

[M, Minv, MiB, MiC, MiD, MiE, F, H, G, L, S] = hdg_compute(master, mesh, app);
 
% compute the index mapping
elcon = mesh.elcon;

% q = getq(MiB, MiC, u, uhat, elcon);
% v = getv(MiD, MiE, u, uhat, elcon);

[AE, FE, Ua, Wa] = assemble(M, MiB, MiC, MiD, MiE, F, H, G, L, S);

uhat = linearsystem(AE, FE, elcon);

u = getu(Wa, Ua, uhat, elcon);
q = getq(MiB, MiC, u, uhat, elcon);
v = getv(MiD, MiE, u, uhat, elcon);

end

function [AE, FE, Ua, Wa] = assemble(M, MiB, MiC, MiD, MiE, F, H, G, L, S)

npv = size(M,1);
[npfe, nt] = size(S);
nd = size(MiB,3);

Qa = zeros(npv*nd,npfe);
Ra = zeros(npv*nd,1);

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
    
    % Di*U = Ei*Uhat + F_K
    Di = -Mass*MiD(:,:,1,1,i);
    Ei = Mass*MiE(:,:,1,1,i);
    for l = 2:nd
      Di = Di - Mass*MiD(:,:,l,l,i);
      Ei = Ei + Mass*MiE(:,:,l,l,i);
    end        
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


