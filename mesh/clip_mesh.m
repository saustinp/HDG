function [p,t] = clip_mesh(p,t,Xc,Nc)
%CLIP_MESH - Clip the mesh, removing all elements on one side of a plane
% That plane is defined by the point Xc and the normal Nc. That normal
% defines a 'positive' side of the plane. Each node that is on the negative
% side is removed from array p. Each element that has at least one node on
% the negative side of the plane is removed from array t.
%
% SYNTAX:  [p,t] = clip_mesh(p,t,Xc,Nc)
%
% INPUTS:
%    p  - Array of nodes coordinates
%    t  - Element to nodes connectivity
%    Xc - (1x3 array) position of one node of the cutting plane
%    Nc - (1x3 array) coordinates of the plane's normal
%
% OUTPUTS:
%    p - Reduced list of nodes coordinates
%    t - Reduced list of elements
%
% Author(s): HDG Team - best team ever.
% Massachusetts Institute of Technology, Dept. of Aeronautics and Astronautics
% email address: XXXXXXX@mit.edu 
% Website: http://aeroastro.mit.edu/
% January 2018; Last revision: January 2018

% Some dimensions
np = size(p,1);
ne = size(t,1);

% Check input format of Xc and Nc
if ~isequal(size(Xc),[1 3])
    error('Error : input 3 has the wrong dimensions.')
end
if ~isequal(size(Nc),[1 3])
    error('Error : input 4 has the wrong dimensions.')
end

% Find suppressed nodes at one side of the plane
MA   = p - ones(size(p,1),1)*Xc;
Dcut = Nc(1)*MA(:,1) + Nc(2)*MA(:,2) + Nc(3)*MA(:,3);
Dcut = Dcut>0;

% Index of nodes suppression for p
j = 1; indp = zeros(np,1);
for i=1:np
    if Dcut(i)
        indp(i) = j;
        j = j+1;
    end
end
% Reduced set of nodes
p = p(Dcut,:);

% Index of element suppression for t
indt = false(ne,1);
for i=1:ne
    indt(i) = all(Dcut(t(i,:)));
end
% Reduced set of elements
t = t(indt,:);

% reindexing the nodes
t = indp(t);

end