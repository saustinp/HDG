function [M, B, F, Minv, P] = hdg_volint(master, app, dgnodes, UDG)

ne   = size(dgnodes,3);
nd   = master.nd;
npv  = master.npv;
ngv  = master.ngv;

% Shap functions and derivatives
shapvt = master.shapvt;
shapvg = master.shapvg(:,:,1);
shapvgdotshapvl = reshape(master.shapvgdotshapvl,[npv*npv ngv (nd+1)]);

% DG solution at Gauss points
nc = size(UDG,2);
udgg = reshape(permute(UDG,[1 3 2]),[npv ne*nc]);
udgg = shapvt(:,:,1)*udgg;
udgg = reshape(udgg,[ngv*ne nc]);

% Positions, Jacobian, and its determinant at Gauss points  
[pg, Xx, jac] = volgeom(master.shapmv, permute(dgnodes,[1 3 2]));

M = reshape(shapvgdotshapvl(:,:,1)*reshape(jac,[ngv ne]),[npv npv ne]);

% convection matrix (u, nabla w)_K
B = zeros(npv,npv,nd,ne);
for i=1:nd
  tmp = reshape(shapvgdotshapvl(:,:,2)*Xx(:,:,i,1),[npv npv ne]);
  for j=2:nd
      tmp = tmp + reshape(shapvgdotshapvl(:,:,j+1)*Xx(:,:,i,j),[npv npv ne]);
  end    
  B(:,:,i,:) = reshape(tmp,[npv, npv, 1, ne]);    
end        

% vector due to source term (s, w)_K
source   = str2func(app.source);
s = source( pg, udgg, app.param); 
F = shapvg*reshape(s.*jac,[ngv ne]);
P = shapvg*reshape(jac,[ngv ne]);

% inverse of the mass matrix
if nargout>3
  Minv = 0*M;
  for i = 1:ne
    Minv(:,:,i) = inv(M(:,:,i));
  end
end







