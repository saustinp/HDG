function F = hdg_sourcevector(master, app, dgnodes, UDG)

ne   = size(dgnodes,3);
nd   = master.nd;
npv  = master.npv;
ngv  = master.ngv;

% Shap functions and derivatives
shapvt = master.shapvt;
shapvg = reshape(master.shapvg,[npv ngv*(nd+1)]);

% DG solution at Gauss points
nc = size(UDG,2);
udgg = reshape(permute(UDG,[1 3 2]),[npv ne*nc]);
udgg = shapvt(:,:,1)*udgg;
udgg = reshape(udgg,[ngv*ne nc]);

% Positions, Jacobian, and its determinant at Gauss points  
[pg, ~, jac] = volgeom(master.shapmv,dgnodes);

% vector due to source term (s, w)_K
s = source( pg, udgg, app.param); 
F = shapvg*reshape(s.*jac,[ngv ne]);

