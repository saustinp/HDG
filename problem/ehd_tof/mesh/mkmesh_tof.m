function mesh = mkmesh_tof(porder, fname)

[p,t] = gmshcall(fname, 2, 0);
r_tip = 5e-4;
p = p/r_tip;    % Nondimensionalization

% % Below is nondimensional
bndexpr = {'all(p(:,1)<min(p0(:,1))+1e-6*5e-4)',...
           'all(p(:,1)>max(p0(:,1))-1e-6*5e-4)',...
           'all(p(:,2)<min(p0(:,2))+1e-6*5e-4)',...
           'all(p(:,2)>max(p0(:,2))-1e-6*5e-4)'};

% Below is dimensional
% bndexpr = {'all(p(:,1)<min(p0(:,1))+1e-9)',...
%            'all(p(:,1)>max(p0(:,1))-1e-9)',...
%            'all(p(:,2)<min(p0(:,2))+1e-9)',...
%            'all(p(:,2)>max(p0(:,2))-1e-9)'};

% Boundaries
% 1. Symmetry boundary
% 2. Outside of cylinder
% 3. Start electrode
% 4. End electrode

elemtype = 0;
nodetype = 1;
mesh = mkmesh(p',t', porder, bndexpr, elemtype, nodetype);
% 
% figure(1); clf;
% meshplot(mesh);
% hold on;
% boundaryplot(mesh,0,1); pause;
% boundaryplot(mesh,0,2); pause;
% boundaryplot(mesh,0,3); pause;
% boundaryplot(mesh,0,4);