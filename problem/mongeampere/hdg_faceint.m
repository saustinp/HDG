function [C, C1, C2, D, E, G, L, H, S] = hdg_faceint(master, mesh, app, dgnodes, UDG)
% FACEINTND compute face integrals 

% Obtain dgnodes, Jacobian matrix and determinant at Gauss points
[pg, nlg, jac] = facegeom(master.shapmf, permute(dgnodes,[1 3 2]), master.permgeom);

ne   = size(dgnodes,3);
nd   = master.nd;
npv  = master.npv;
npf  = master.npf;
ngf  = master.ngf;
nfe  = size(master.perm,2);
npfe = npf*nfe;

% Shap functions 
perm            = master.perm(:,:,1);
shapft          = master.shapft(:,:,1);
shapfg          = master.shapfg(:,:,1);
shapfgdotshapfc = master.shapfgdotshapfc(:,:,1);

C = zeros(npv,npfe,nd,ne);  % <uhat, w n>_{\partial K}              
G = zeros(npfe,npv,nd,ne);  % <q dot n, mu>_{\partial K}                    
for i=1:nd    
    njc = reshape(nlg(:,i).*(jac),[ngf,nfe*ne]);
    wrk = reshape(shapfgdotshapfc*njc,[npf npfe 1 ne]);   
    for j=1:nfe
        C(perm(:,j),(j-1)*npf+1:j*npf,i,:) = C(perm(:,j),(j-1)*npf+1:j*npf,i,:)  + wrk(1:npf,(j-1)*npf+1:j*npf,1,:);      
        G((j-1)*npf+1:j*npf,perm(:,j),i,:) = G((j-1)*npf+1:j*npf,perm(:,j),i,:)  + wrk(1:npf,(j-1)*npf+1:j*npf,1,:);      
    end    
end

C1 = zeros(npv,npv,nd,ne);   % <u, w n*n>_{\partial K}                   
C2 = zeros(npv,npfe,nd,ne);   % <uhat, w n*n>_{\partial K}                   
for i=1:nd    
    njc = reshape(nlg(:,i).*nlg(:,i).*(jac),[ngf,nfe*ne]);
    wrk = reshape(shapfgdotshapfc*njc,[npf npfe 1 ne]);   
    for j=1:nfe
        C1(perm(:,j),perm(:,j),i,:) = wrk(1:npf,(j-1)*npf+1:j*npf,1,:);      
        C2(perm(:,j),(j-1)*npf+1:j*npf,i,:) = wrk(1:npf,(j-1)*npf+1:j*npf,1,:);      
    end          
end

% mass matrices on faces
njc = reshape(jac,[ngf,nfe*ne]);
wrk = reshape(shapfgdotshapfc*njc,[npf npfe ne]);  

% matrices from mass matrices on faces
D = zeros(npv,npv,ne);   % <u, w>_{\partial K}        
E = zeros(npv,npfe,ne);  % <uhat, w>_{\partial K}              
L = zeros(npfe,npv,ne);  % <u, \mu>_{\partial K}              
H = zeros(npfe,npfe,ne); % <eta, \mu>_{\partial K}                     
for j=1:nfe
    D(perm(:,j),perm(:,j),:) = D(perm(:,j),perm(:,j),:) + wrk(1:npf,(j-1)*npf+1:j*npf,:);      
    E(perm(:,j),(j-1)*npf+1:j*npf,:) = E(perm(:,j),(j-1)*npf+1:j*npf,:) + wrk(1:npf,(j-1)*npf+1:j*npf,:);      
    L((j-1)*npf+1:j*npf,perm(:,j),:) = L((j-1)*npf+1:j*npf,perm(:,j),:) + wrk(1:npf,(j-1)*npf+1:j*npf,:);      
    H((j-1)*npf+1:j*npf,(j-1)*npf+1:j*npf,:) = H((j-1)*npf+1:j*npf,(j-1)*npf+1:j*npf,:)  + wrk(1:npf,(j-1)*npf+1:j*npf,:);      
end

nlg = reshape(nlg, [ngf nfe ne nd]);
pg = reshape(pg, [ngf nfe ne nd]);
jac = reshape(jac, [ngf nfe ne]);

% vector due to the boundary data
fbou   = str2func(app.fbou);
S = zeros(npfe,ne);  % <g, \mu>_{\partial K}              

for i = 1:ne % loop over each element
  for j = 1:nfe % loop over each face of the element i
    ibf = mesh.f(abs(mesh.t2f(i,j)),end);
    if (ibf<0) % boundary face
      ibc = app.bcm(-ibf);    
      I = perm(:,j,1);
      J = ((j-1)*npf+1):j*npf;            
            
      ugf = shapft*UDG(I,:,i);
      nlf = reshape(nlg(:,j,i,:),[ngf nd]);
      xgf = reshape(pg(:,j,i,:),[ngf nd]);                          
                        
      % vector due to boundary data
      gb = fbou(xgf, nlf, ugf, app.param, ibc);                                             
      S(J,i) = S(J,i) + shapfg*(jac(:,j,i).*gb); 
      
      if (ibc == 0) % Dirichlet boundary condition       
        G(J,I,:,i) = 0;
        L(J,I,i) = 0;                
      elseif (ibc > 1) % Nonstandard Neumann boundary condition       
        G(J,I,:,i) = 0;
      end
    end
  end
end



