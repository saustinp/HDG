function [Ae, Re, De, We, Ue, MiB, MiC, MiD, MiE, G, L, H, S, P] = hdg_precompute(master, mesh, app, UDG)

[M, B, F, Minv, P] = hdg_volint(master, app, mesh.dgnodes, UDG);

[C, C1, C2, D, E, G, L, H, S] = hdg_faceint(master, mesh, app, mesh.dgnodes,UDG);

ne = size(mesh.dgnodes,3);
nd   = master.nd;
npv  = master.npv;
npf  = master.npf;
nfe  = size(master.perm,2);
npfe = npf*nfe;

MiB = zeros(npv,npv,nd,ne);
MiC = zeros(npv,npfe,nd,ne);
MiD = zeros(npv,npv,nd,nd,ne);
MiE = zeros(npv,npfe,nd,nd,ne);
We = zeros(npv, ne);
Ue = zeros(npv, npfe, ne);
Re = zeros(npfe, ne);
Ae = zeros(npfe, npfe, ne);
De = zeros(npv, npv, ne);

MiBt = zeros(npv,npv,nd);
for i = 1:ne
    Mi = Minv(:,:,i);        
    % q - nabla u = 0 
    % -> (q, w) + (u, nabla dot w) - <uhat, w dot n> = 0 
    % -> M q + B*u - C*uhat = 0 
    % -> q = -(Minv*B)*u + (Minv*C)*uhat  
    for l = 1:nd      
      MiB(:,:,l,i) = Mi*B(:,:,l,i);
      MiBt(:,:,l) = Mi*B(:,:,l,i)';
      MiC(:,:,l,i) = Mi*C(:,:,l,i);
    end            
    
    % v - nabla q = 0    
    % -> (v, w) + (q, nabla w) - <qhat, w dot n> = 0 
    % -> M v + B*q - C*qhat = 0
    % -> (v, w) - (nabla q, w) - <(u - uhat) tau n, w dot n> = 0     
    % -> M v - B'*q + C1*u - C2*uhat = 0
    % -> v = (Minv*B')*q - (Minv*C1)*u + (Minv*C2)*uhat 
    % -> v = (Minv*B')*(-(Minv*B)*u + (Minv*C)*uhat) - (Minv*C1)*u + (Minv*C2)*uhat 
    % v = -((Minv*B')*(Minv*B) + (Minv*C1))*u + ((Minv*B')*(Minv*C) + (Minv*C2))*uhat    
    if nd==2
        MiD(:,:,1,1,i) = -MiBt(:,:,1)*MiB(:,:,1,i) - Mi*C1(:,:,1,i);
        MiD(:,:,2,2,i) = -MiBt(:,:,2)*MiB(:,:,2,i) - Mi*C1(:,:,2,i);
        MiD(:,:,1,2,i) = -MiBt(:,:,2)*MiB(:,:,1,i) + Mi*C1(:,:,2,i);
        MiD(:,:,2,1,i) = -MiBt(:,:,1)*MiB(:,:,2,i) + Mi*C1(:,:,1,i);
        
        MiE(:,:,1,1,i) = MiBt(:,:,1)*MiC(:,:,1,i) + Mi*C2(:,:,1,i);      
        MiE(:,:,2,2,i) = MiBt(:,:,2)*MiC(:,:,2,i) + Mi*C2(:,:,2,i);      
        MiE(:,:,1,2,i) = MiBt(:,:,2)*MiC(:,:,1,i) - Mi*C2(:,:,2,i);      
        MiE(:,:,2,1,i) = MiBt(:,:,1)*MiC(:,:,2,i) - Mi*C2(:,:,1,i);      
    end
    
    % - nabla dot q = s
    % (q, nabla w) - <qhat, w n>  = (s, w)
    %  qhat = q - tau (u - uhat) n
    % -(nabla dot q,  w) + <tau u, w> - <tau uhat, w>   = (s, w)
    % -B'*q + D*u - E*uhat = F
    % -B'*(-(Minv*B)*u + (Minv*C)*uhat) + D*u - E*uhat = F
    % (B'*Minv*B + D)*u - (B'*Minv*C + E) *uhat = F
    % u =  inv(B'*Minv*B + D)*(F + (B'*Minv*C + E)*uhat) 
    % We = inv(B'*Minv*B + D)*F
    % Ue = inv(B'*Minv*B + D)*(B'*Minv*C + E)
    % u =  We + Ue*uhat 
    
    % Compute B'*Minv*B + D and B'*Minv*C + E
    De(:,:,i) = D(:,:,i);
    Ee = E(:,:,i);
    for l = 1:nd      
      De(:,:,i) = De(:,:,i) + B(:,:,l,i)'*MiB(:,:,l,i);      
      Ee = Ee + B(:,:,l,i)'*MiC(:,:,l,i);      
    end            
        
    % compute We and Ue
    We(:,i) = De(:,:,i)\F(:,i);  
    Ue(:,:,i) = De(:,:,i)\Ee;    
     
    % <qhat dot n, mu> + <f(q, u, uhat, n, UDG), mu> = <g, mu> 
    % -> G*q - L*u + H*uhat = S
    % -> G*(-(Minv*B)*u + (Minv*C)*uhat)  - L*u + H*uhat = S
    % -> G*(-(Minv*B)*(We + Ue*uhat)  + (Minv*C)*uhat)  - L*(We + Ue*uhat)  + H*uhat = S
    % -> G*(Minv*C - (Minv*B)*Ue)*uhat - L*Ue*uhat  + H*uhat = S + G*(Minv*B)*We + L*We
    % -> Re = S + G*(Minv*B)*We + L*We = S + (L + G*(Minv*B))*We
    % -> Ae = G*(Minv*C - (Minv*B)*Ue) - L*Ue + H = H + G*Minv*C - (L + G*(Minv*B))*Ue     
    
    % G = zeros(npfe,npv,nd,ne);  % <q dot n, mu>_{\partial K}                    
    % Compute G*(Minv*B) and G*(Minv*C)
    GMiB = G(:,:,1,i)*MiB(:,:,1,i);    
    GMiC = G(:,:,1,i)*MiC(:,:,1,i);    
    for l = 2:nd
      GMiB = GMiB + G(:,:,l,i)*MiB(:,:,l,i);    
      GMiC = GMiC + G(:,:,l,i)*MiC(:,:,l,i);    
    end
    
    Re(:,i) = S(:,i) + (L(:,:,i) + GMiB)*We(:,i);       
    Ae(:,:,i) = H(:,:,i) + GMiC - (L(:,:,i) + GMiB)*Ue(:,:,i);                
end


