function [ dmin ] = get_min_dist( dgnodes )
%GET_MIN_DIST - compute the lowest distance between two DG-mesh nodes
%This function is usefull for the computation of the global CFL number.
%
% SYNTAX:  [dmin] = get_min_dist( dgnodes )
%
% INPUTS:
%    dgnodes - 3D array of DG-nodes position ( it is mesh.dgnodes field)
%
% OUTPUTS:
%    dmin - minimal distance between nodes
%
% OTHER M-FILES REQUIRED: none
% SUBFUNCTIONS: none
% MAT-FILES REQUIRED: none
%
% SEE ALSO: COMPUTE_CFL

% Author(s): Sebastien Terrana
% Massachusetts Institute of Technology, Dept. of Aeronautics and Astronautics
% email address: terrana@mit.edu 
% Website: http://aeroastro.mit.edu/
% January 2018; Last revision: January 2018

%------------- BEGIN CODE --------------
% brute-force search on all element nodes, for each element
[npv,nd,ne] = size(dgnodes);
dminEl = realmax*ones(ne,1);
for i=1:npv
    for j=i+1:npv
        distmp = squeeze(sqrt((dgnodes(j,1,:)-dgnodes(i,1,:)).^2 ...
                             +(dgnodes(j,2,:)-dgnodes(i,2,:)).^2 ...
                             +(dgnodes(j,3,:)-dgnodes(i,3,:)).^2));
        dminEl = min(dminEl,distmp);
    end
end
% Return the global minimum distance
dmin = min(dminEl);
end