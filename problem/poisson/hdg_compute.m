function [M, Minv, MiB, MiC, MiD, MiE, F, H, G, L, S] = hdg_compute(master, mesh, app)

tau = app.tau;
param = app.param;
bcm = app.bcm;
fbou   = str2func(app.fbou);
source   = str2func(app.source);

nd = mesh.nd;
nt  = size(mesh.t,1);
npv = size(mesh.dgnodes,1);
npf = size(master.plocfc,1);
nfe = size(master.perm,2);
ngf = master.ngf;

shap = squeeze(master.shapvl(:,:,1));
% shapxi = squeeze(master.shapvl(:,:,2))*diag(master.gwvl);
% shapet = squeeze(master.shapvl(:,:,3))*diag(master.gwvl); 
shapfc = squeeze(master.shapfc(:,:,1));
perm = master.perm;

% allocate memory for the elemental matrices
B_K = zeros(nd*npv,npv);
MiBt = zeros(npv,npv,nd);
tmn = zeros(npf,npf,nd);
tm2 = zeros(npf,npf,nd);

M = zeros(npv,npv,nt);
Minv = zeros(npv,npv,nt);
MiB = zeros(npv,npv,nd,nt);
MiC = zeros(npv,npf*nfe,nd,nt);
MiD = zeros(npv,npv,nd,nd,nt);
MiE = zeros(npv,npf*nfe,nd,nd,nt);
F = zeros(npv,nt);
H = zeros(nfe*npf,nfe*npf,nt);    
G = zeros(nfe*npf,nd*npv,nt);
L = zeros(nfe*npf,npv,nt);
S = zeros(nfe*npf,nt); 

for i = 1:nt % loop over each element
    [xg, Xx, jac] = volgeom(master.shapmv,reshape(mesh.dgnodes(:,:,i),[npv 1 nd]));      
    
    Mass  = shap*diag(master.gwvl.*jac)*shap';
    M(:,:,i) = Mass;
    Minv(:,:,i) = inv(Mass);
        
    for l=1:nd
        tmp = reshape(master.shapvgdotshapvl(:,:,2)*Xx(:,:,l,1),[npv npv]);
        for j=2:nd
            tmp = tmp + reshape(master.shapvgdotshapvl(:,:,j+1)*Xx(:,:,l,j),[npv npv]);
        end
        B_K(((l-1)*npv+1):l*npv,1:npv) = tmp;
    end                        
    
    % form F_K
    fg = source(xg, [xg xg(:,1)], param);
    F(:,i) = shap*(master.gwvl.*jac.*fg);
            
    % allocate memory for the elemental matrices
    D1_K = zeros(npv,npv,nd);
    E1_K = zeros(npv,nfe*npf,nd);
    C_K = zeros(nd*npv,nfe*npf);
    D_K = zeros(npv,npv);
    E_K = zeros(npv,nfe*npf);
    M_K = zeros(nfe*npf,nfe*npf);    
    G_K = zeros(nfe*npf,nd*npv);
    I_K = zeros(nfe*npf,npv);
    S_K = zeros(nfe*npf,1); 
    
    [xf, nlg, jacf] = facegeom(master.shapmf,reshape(mesh.dgnodes(:,:,i),[npv 1 nd]),perm);        
    xf = reshape(xf, [ngf nfe nd]);
    nlg = reshape(nlg, [ngf nfe nd]);
    jacf = reshape(jacf, [ngf nfe]);
    
    for j = 1:nfe % loop over each face of the element i
        I = perm(:,j,1);
        J = ((j-1)*npf+1):j*npf;        
        
        nl = reshape(nlg(:,j,:),[ngf nd]);
        xgf = reshape(xf(:,j,:),[ngf nd]);
        dws = master.gwfc.*jacf(:,j);                        
        
        % form C_K
        for l=1:nd
            tmn(:,:,l) = shapfc*diag(dws.*nl(:,l))*shapfc';
            C_K((l-1)*npv+I,J) = C_K((l-1)*npv+I,J) + tmn(:,:,l);   
            tm2(:,:,l) = shapfc*diag(dws.*nl(:,l).*nl(:,l))*shapfc';
            D1_K(I,I,l) = D1_K(I,I,l) + tau*tm2(:,:,l);   
            E1_K(I,J,l) = E1_K(I,J,l) + tau*tm2(:,:,l);   
        end

        tmp = shapfc*diag(tau*dws)*shapfc';
        
        % form D_K
        D_K(I,I) = D_K(I,I) + tmp;

        % form E_K
        E_K(I,J) = E_K(I,J) + tmp;

        % form M_K
        M_K(J,J) = M_K(J,J) + tmp;        
               
        ibf = mesh.f(abs(mesh.t2f(i,j)),end);
        if ibf>0            
            % form G_K
            for l=1:nd
                G_K(J,(l-1)*npv+I) = G_K(J,(l-1)*npv+I) + tmn(:,:,l);   
            end
            
            % form I_K
            I_K(J,I) = I_K(J,I) + tmp;                        
        else
            ibc = bcm(-ibf);
            gb = fbou(xgf,nl,param,-ibf);
            S_K(J,1) = S_K(J,1) + shapfc*(dws.*gb); 
                        
            if (ibc==1) % Neumman
                % form G_K
                for l=1:nd
                    G_K(J,(l-1)*npv+I) = G_K(J,(l-1)*npv+I) + tmn(:,:,l);   
                end

                % form I_K
                I_K(J,I) = I_K(J,I) + tmp;                                        
            end            
        end        
    end    
        
    H(:,:,i) = M_K;    
    G(:,:,i) = G_K;
    L(:,:,i) = I_K;
    S(:,i) = S_K;
    
    % compute: Minv, B, C, D, E, F, H, G, L, S
    % Store Minv, MiB, MiC, MiD, MiE, F, H, G, L, S

    % q - du/dx = 0 -> M*Q + B*U - C*Uhat = 0
    % v - dq/dx = 0 -> M*V - B'*Q + D1*U - E1*Uhat = 0 
    % V = Mi*B'*(Mi*C*uhat - Mi*B*U)) - Mi*D1*U + Mi*E*Uhat
    % V = -(Mi*B'*Mi*B + Mi*D1)*U + (Mi*B'*Mi*C + Mi*E1)*Uhat        
    % V = MiD*U + MiE*Uhat
    
    Mi = Minv(:,:,i);        
    for l = 1:nd      
      MiB(:,:,l,i) = Mi*B_K(((l-1)*npv+1):l*npv,:);
      MiBt(:,:,l) = Mi*B_K(((l-1)*npv+1):l*npv,:)';
      MiC(:,:,l,i) = Mi*C_K(((l-1)*npv+1):l*npv,:);      
    end            
    if nd==2
        MiD(:,:,1,1,i) = -MiBt(:,:,1)*MiB(:,:,1,i) - Mi*D1_K(:,:,1);
        MiD(:,:,2,2,i) = -MiBt(:,:,2)*MiB(:,:,2,i) - Mi*D1_K(:,:,2);
        MiD(:,:,1,2,i) = -MiBt(:,:,2)*MiB(:,:,1,i) - Mi*D1_K(:,:,2);
        MiD(:,:,2,1,i) = -MiBt(:,:,1)*MiB(:,:,2,i) - Mi*D1_K(:,:,1);
        
        MiE(:,:,1,1,i) = Mi*E1_K(:,:,1) + MiBt(:,:,1)*MiC(:,:,1,i);      
        MiE(:,:,2,2,i) = Mi*E1_K(:,:,2) + MiBt(:,:,2)*MiC(:,:,2,i);      
        MiE(:,:,1,2,i) = Mi*E1_K(:,:,2) + MiBt(:,:,2)*MiC(:,:,1,i);      
        MiE(:,:,2,1,i) = Mi*E1_K(:,:,1) + MiBt(:,:,1)*MiC(:,:,2,i);      
    end
%     for l = 1:nd      
%       for m = 1:nd      
%         MiD(:,:,l,m,i) = -MiBt(:,:,l)*MiB(:,:,m,i) - (1/nd)*Mi*D_K;
%         MiE(:,:,l,m,i) = (1/nd)*Mi*E_K + MiBt(:,:,l)*MiC(:,:,m,i);
%       end
%     end        
end


