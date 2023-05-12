function [S, G] = hdg_boundaryvector(master,app,dgnodes,UDG,G)
% FACEINTND compute face integrals 

ne   = size(dgnodes,3);
nd   = master.nd;
npf  = master.npf;
nfe  = size(master.perm,2);
npfe = npf*nfe;

% Shap functions 
perm            = master.perm(:,:,1);
shapft          = master.shapft(:,:,1);
shapfg          = master.shapfg(:,:,1);
shapfc          = master.shapfc(:,:,1);

% vector due to the boundary data
fbou   = str2func(app.fbou);
dfbou   = str2func(app.dfbou);
S = zeros(npfe,ne);  % <g, \mu>_{\partial K}              

for i = 1:ne % loop over each element
  for j = 1:nfe % loop over each face of the element i
    ibf = mesh.f(abs(mesh.t2f(i,j)),end);
    if (ibf<0) % boundary face
      ibc = bcm(-ibf);    
      I = perm(:,j,1);
      J = ((j-1)*npf+1):j*npf;            
            
      % Obtain dgnodes, Jacobian matrix and determinant at Gauss points
      ugf = shapft*UDG(I,:,i);
      [xgf, nlf, jacf] = facegeom(master.shapmf,reshape(dgnodes(I,:,i),[npf 1 nd]),master.perm);        
                        
      % vector due to boundary data
      gb = fbou(xgf, ugf, nlf, app.param, ibc);                                             
      S(J,i) = S(J,i) + shapfg*(jacf.*gb);       
      
      if narout>1
        qb = dfbou(xgf, ugf, nlf, app.bcd(-ibf,:), ibc);      
        G(J,:,i) = 0;
        for l=1:nd              
            tmx = shapfc*diag(master.gwfc.*jacf.*qb(:,l))*shapfc';
            G(J,(l-1)*npv+I,i) = G(J,(l-1)*npv+I,i) + tmx;   
        end                        
      end
    end
  end
end




