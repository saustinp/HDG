function mesh = mkmesh_chen(porder, fname)

[p,t] = gmshcall(fname, 2, 0);
r_tip = 220e-6;
p = p/r_tip;    % Nondimensionalization

% Below is nondimensional
bndexpr = {'all(p(:,1)<min(p0(:,1))+1e-6*220e-6)',...
           'all(p(:,1)>max(p0(:,1))-1e-6*220e-6)', ...
           'all(p(:,2)>max(p0(:,2))-1e-6*220e-6)',...
           'all(p(:,2)<min(p0(:,2))+1e-6*220e-6) & all(p(:,1)<0.015*220e-6+1e-6*220e-6)',...
           'all(p(:,2)<min(p0(:,2))+1e-6*220e-6) & all(p(:,1)>0.017*220e-6-1e-6*220e-6)',...
           'all(p(:,2)<-0.01*220e-6)',...
           'all(p(:,2)>-0.01*220e-6)'}; 

% Below is dimensional
% bndexpr = {'all(p(:,1)<min(p0(:,1))+1e-8)',...
%            'all(p(:,1)>max(p0(:,1))-1e-8)', ...
%            'all(p(:,2)>max(p0(:,2))-1e-8)',...
%            'all(p(:,2)<min(p0(:,2))+1e-8) & all(p(:,1)<0.015+1e-8)',...
%            'all(p(:,2)<min(p0(:,2))+1e-8) & all(p(:,1)>0.017-1e-8)',...
%            'all(p(:,2)<-0.01)',...
%            'all(p(:,2)>-0.01)'}; 

% Boundaries
% 1. Axisymmetry
% 2. Right farfield
% 3. Top farfield
% 4. Outflow
% 5. Ground plane
% 6. Cylinder
% 7. Needle
         
elemtype = 0;
nodetype = 1;
mesh = mkmesh(p',t', porder, bndexpr, elemtype, nodetype);

% figure(1); clf;
% meshplot(mesh);
% hold on;
% boundaryplot(mesh,0,1); 
% boundaryplot(mesh,0,2);
% boundaryplot(mesh,0,3);
% boundaryplot(mesh,0,4);
% boundaryplot(mesh,0,5);
% boundaryplot(mesh,0,6);
% boundaryplot(mesh,0,7);





