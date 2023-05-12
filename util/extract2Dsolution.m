
function [mesh2d, UDG2d, additionalFields2d] = extract2Dsolution(mesh,UDG,additionalFields,zTarget,getOnly2Dfields,plotFlag)

% Description: Extract 2D mesh and solution at zTarget from 3D mesh and
% solution. Only valid for quasi-2D problem with extrusion in Z directon.

if nargin < 3; additionalFields = []; end
if nargin < 4; zTarget = 0; end
if nargin < 5; getOnly2Dfields = 1; end
if nargin < 6; plotFlag = 0; end
if mesh.nd ~= 3; error('Input mesh and UDG must be three-dimensional.'); end

face2p_tet = [2,3,4;
              1,4,3;
              1,2,4;
              1,3,2];

face2p_hex = [1,4,3,2;
              5,6,7,8;
              1,2,6,5;
              3,4,8,7;
              2,3,7,6;
              4,1,5,8];

tol = 1e-6;
elemtype = mesh.elemtype;

pInPlane = find((mesh.p(:,3) < zTarget+tol).*(mesh.p(:,3) > zTarget-tol));

% Determine if target plane is the end (in z direction) of computational domain:
if abs(zTarget - max(mesh.p(:,3))) < tol
    lastPlane = 1;
else
    lastPlane = 0;
end

if elemtype == 0
    elemPlane = find(sum(ismember(mesh.t,pInPlane(:)),2) >= 3);
elseif elemtype == 1
    elemPlane = find(sum(ismember(mesh.t,pInPlane(:)),2) >= 4);
end

elemPlane_New = [];
if elemtype == 0; mesh2d.t = zeros(1,3);
elseif elemtype == 1; mesh2d.t = zeros(1,4); end
UDG2d = zeros(size(mesh.perm,1),size(UDG,2),1);
dgnodes2d = zeros(size(mesh.perm,1),2,1);
if ~isempty(additionalFields); additionalFields2d = zeros(size(mesh.perm,1),size(additionalFields,2),1); end
zCoord = mesh.p(:,3);
for elem = elemPlane(:)'
    psInElem = mesh.t(elem,:);
    if (~lastPlane && all(mesh.p(psInElem(:),3) > zTarget - tol)) || (lastPlane && all(mesh.p(psInElem(:),3) < zTarget + tol))
        elemPlane_New = [elemPlane_New; elem];
        if lastPlane
            localFace = find(all(zCoord(mesh.f(mesh.t2f(elem,:),1:end-2)) > zTarget-tol, 2));
        else
            localFace = find(all(zCoord(mesh.f(mesh.t2f(elem,:),1:end-2)) < zTarget+tol, 2));
        end
        if length(localFace) ~= 1; error('Something wrong'); end
        if elemtype == 0; mesh2d.t(length(elemPlane_New),:) = mesh.t(elem,face2p_tet(localFace,:));
        elseif elemtype == 1; mesh2d.t(length(elemPlane_New),:) = mesh.t(elem,face2p_hex(localFace,:)); end
%         if elemtype == 0; mesh2d.t(length(elemPlane_New),:) = mesh.f(mesh.t2f(elem,localFace),[1,2,3]);
%         elseif elemtype == 1; mesh2d.t(length(elemPlane_New),:) = mesh.f(mesh.t2f(elem,localFace),[1,2,3,4]); end
        UDG2d(:,:,length(elemPlane_New)) = UDG(mesh.perm(:,localFace),:,elem);
        dgnodes2d(:,:,length(elemPlane_New)) = mesh.dgnodes(mesh.perm(:,localFace),1:2,elem);
        if ~isempty(additionalFields); additionalFields2d(:,:,length(elemPlane_New)) = additionalFields(mesh.perm(:,localFace),:,elem); end
    end
end
elemPlane = elemPlane_New;

% Define 2d vertex locations:
verticesInPlane = sort(unique(mesh2d.t(:)));
mesh2d.p = mesh.p(verticesInPlane,1:2);

% Renumber the 2d vertices from 1 to mesh2d.np:
map = zeros(1,max(mesh2d.t(:)));
mesh2d_t_unique = unique(mesh2d.t);
map(sort(mesh2d_t_unique)) = 1:length(mesh2d_t_unique);
mesh2d.t = map(mesh2d.t);
if min(mesh2d.t(:)) == 0; error('Something wrong'); end

[mesh2d.p,mesh2d.t]=fixmesh(mesh2d.p,mesh2d.t);
mesh2d = mkmesh(mesh2d.p,mesh2d.t,mesh.porder,{'true'},elemtype,0);
mesh2d.dgnodes = dgnodes2d;

if getOnly2Dfields == 1; UDG2d = UDG2d(:,[1,2,3,5,6,7,8,9,11,12,13,15],:); end

if plotFlag; scaplot(mesh2d,UDG2d(:,2,:)); end

end
