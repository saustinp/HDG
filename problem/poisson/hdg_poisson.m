function [u,q,uhat] = hdg_poisson(master, mesh, app)
%HDG_POI Solve the poission problem using HDG
%
%   [u,q,uhat] = HDG_POI(MASTER,MESH,kappa, tau)
%
%      MASTER:       Master structure
%      MESH:         Mesh structure
%      kappa:        Diffusion coefficient
%      tau:          Stabilization parameter
%

kappa = app.kappa;
tau = app.tau;
coord = app.coord;
param = app.param;
bcm = app.bcm;
fbou   = str2func(app.fbou);
source   = str2func(app.source);

ne = mesh.ne;
nd = mesh.nd;
nt  = size(mesh.t,1);
%nf  = size(mesh.f,1);
npv = size(mesh.dgnodes,1);
npf = size(master.plocfc,1);
nfe = size(master.perm,2);
ngf = master.ngf;

shap = squeeze(master.shapvl(:,:,1));
% shapxi = squeeze(master.shapvl(:,:,2))*diag(master.gwvl);
% shapet = squeeze(master.shapvl(:,:,3))*diag(master.gwvl); 
shapfc = squeeze(master.shapfc(:,:,1));
perm = master.perm;

% compute the index mapping
%elcon = elconnectivities(npf,mesh.t2f);
elcon = mesh.elcon;

% allocate memory for the elemental matrices
A_K = zeros(nd*npv,nd*npv);
B_K = zeros(nd*npv,npv);
tmn = zeros(npf,npf,nd);

% allocate memory for the global matrix and vector
% H = sparse(nf*npf,nf*npf);
% R = sparse(nf*npf,1);

% allocate memory 
% P = zeros((nd+1)*npv,nfe*npf,nt);
% L = zeros((nd+1)*npv,nt);
U = zeros(npv,nfe*npf,nt);
W = zeros(npv,nt);
Minv = zeros(npv,npv,nt);
B = zeros(npv*nd,npv,nt);
C = zeros(npv*nd,npf*nfe,nt);

AE = zeros(npf*nfe,npf*nfe,nt);
FE = zeros(npf*nfe,nt);
for i = 1:nt % loop over each element
    [xg, Xx, jac] = volgeom(master.shapmv,reshape(mesh.dgnodes(:,:,i),[npv 1 nd]));      
    
    if coord==0
        Mass  = shap*diag(master.gwvl.*jac)*shap';
    else
        Mass  = shap*diag(master.gwvl.*jac./xg(:,1))*shap';   
    end
    Minv(:,:,i) = inv(shap*diag(master.gwvl.*jac)*shap');
    
    % form A_K
    for l=1:nd
        A_K(((l-1)*npv+1):l*npv,((l-1)*npv+1):l*npv) = (1/kappa)*Mass;
    end
    
    for l=1:nd
        tmp = reshape(master.shapvgdotshapvl(:,:,2)*Xx(:,:,l,1),[npv npv]);
        for j=2:nd
            tmp = tmp + reshape(master.shapvgdotshapvl(:,:,j+1)*Xx(:,:,l,j),[npv npv]);
        end
        B_K(((l-1)*npv+1):l*npv,1:npv) = tmp;
    end                        
    B(:,:,i) = B_K;    
    
    % form F_K
    %xg = squeeze(master.shapvl(:,:,1))'*mesh.dgnodes(:,:,i);
    fg = source(xg, [xg xg(:,1)], param);
    F_K = shap*(master.gwvl.*jac.*fg);
            
    % allocate memory for the elemental matrices
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
            % form C_K
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
                % form C_K
                for l=1:nd
                    G_K(J,(l-1)*npv+I) = G_K(J,(l-1)*npv+I) + tmn(:,:,l);   
                end

                % form I_K
                I_K(J,I) = I_K(J,I) + tmp;                                        
            end            
        end        
    end    
    C(:,:,i) = C_K;
    
    invA = inv(A_K);        
    U(:,:,i) = ((B_K'*invA)*B_K+D_K)\((B_K'*invA)*C_K+E_K);
    W(:,i) = ((B_K'*invA)*B_K+D_K)\F_K;    
    P = [invA*(C_K-B_K*U(:,:,i)); U(:,:,i)];   
    L = [-invA*B_K*W(:,i); W(:,i)];   

    % these are needed to compute u and q
%     P(:,:,i) = [A_K B_K; -B_K' D_K]\[C_K; E_K];
%     L(:,i) = [A_K B_K; -B_K' D_K]\[zeros(nd*npv,1); F_K];
   
    % form elemental stiffness matrix 
    H_K = M_K + [G_K -I_K]*P;
                
    % form elemental residual vector 
    R_K = -[G_K -I_K]*L; 
    
    % assemble the elemental matrices and vectors into the global matrix and vector
%     imap = elcon(:,i);
%     H(imap,imap) = H(imap,imap) + H_K;
%     R(imap) = R(imap) + R_K + S_K;
    
    AE(:,:,i) = H_K;
    FE(:,i) = R_K + S_K;    
    
    if rem(i,100)==0
        i
    end
end

% if app.directsolver==0
%     disp('Assembly...');
%     nch = 1;
%     [AE,FE] = mkDenseBlockSystem(reshape(AE,[nch npf*nfe nch npf*nfe nt]), reshape(FE,[nch npf*nfe nt]), mesh.f, mesh.t2f, mesh.elcon);
%     f2f = mkf2f(mesh.f, mesh.t2f);
% 
%     y = mesh.dgnodes(:,2,:);
%     ymin = min(y(:));
%     ymax = max(y(:));
%     udg = param(1) - ((y-ymin)/(ymax-ymin))*param(1);
%     uh = inituhat(master,mesh.elcon,udg,1);
%     disp('GMRES...');
%     uhat = gmres(AE, FE, f2f, uh);
%     uhat = uhat(:);
% else        
    % solve for uhat
    il = zeros(nfe*npf,nfe*npf,ne);
    jl = zeros(nfe*npf,nfe*npf,ne);
    for i=1:ne
        con = repmat((elcon(:,i)'-1),1,1)+ones(1,nfe*npf);
        con = reshape(con,nfe*npf,1);
        il(:,:,i) = repmat(con ,1,nfe*npf);
        jl(:,:,i) = repmat(con',nfe*npf,1);        
    end
    H = sparse(reshape(il,(npf*nfe)^2*ne,1),reshape(jl,(npf*nfe)^2*ne,1),reshape(AE,(npf*nfe)^2*ne,1));
    R = sparse(reshape(il(:,1,:),(npf*nfe)*ne,1),ones((npf*nfe)*ne,1),reshape(FE,(npf*nfe)*ne,1));                    
    uhat = H\R;
%end

% compute u and q
u = zeros(npv,1,nt);
q = zeros(npv,2,nt);
for i = 1:nt
    imap = elcon(:,i);
%     QU = L(:,i) + P(:,:,i)*uhat(imap);
%     q(:,:,i) = reshape(QU(1:2*npv),[npv 2]);
%     u(:,1,i) = QU(2*npv+1:end);    
    
    u(:,1,i) = W(:,i) + U(:,:,i)*uhat(imap);
    tmp = -B(:,:,i)*u(:,1,i) + C(:,:,i)*uhat(imap);
    q(:,:,i) = Minv(:,:,i)*reshape(tmp,[npv nd]);    
end

