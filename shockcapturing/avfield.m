function a = avfield(div, cgelcon, cgent2dgent, rowent2elem, nodeR, weightR, alpha)

% velocity divergence
% div = hp.*(UDG(:,2,:) + UDG(:,3,:));

% cg field 
scg = udg2ucg(div, cgent2dgent, rowent2elem);

% filtering cg field
s = filtering(scg, nodeR, weightR);

% dg field
s = ucg2udg(s, cgelcon);

% artificial viscosity field
a = s.*(atan(alpha*s)/pi + 0.5) - atan(alpha)/pi + 0.5;

