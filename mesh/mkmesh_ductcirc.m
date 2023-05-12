function mesh = mkmesh_ductcirc(mesh,rb,H,L)
%MKMESH_DUCT Map a unit square mesh to a cos^2 Duct
%   MESH = MKMESH_DUCT(MESH,DB,DT,H)
%
%      MESH:     Mesh data structure
%                   input: mesh for the unit square created with
%                          mkmesh_square
%                   output: mesh for duct
%      DB:       Height of bottom bump
%      DT:       Height of top bump
%      H:        Height of channel

x1 = (L-2)/2;
x2 = x1+2;

a = 1;
p = mesh.p;
p(:,2) = loginc(p(:,2),a);
pnew(:,1) = L*p(:,1);
pnew(:,2) = H*p(:,2).*(pnew(:,1)<=x1 | pnew(:,1)>=x2) + ...
            (p(:,2).*H + ...
             (1-p(:,2)).*(sqrt(abs(rb^2-(pnew(:,1)-L/2).^2)))).* ...
            (pnew(:,1)>x1 & pnew(:,1)<x2);
mesh.p = pnew;

if isfield(mesh,'dgnodes') && ~isempty(mesh.dgnodes)
   clear pnew;
   p = mesh.dgnodes;
   p(:,2,:) = loginc(p(:,2,:),a);
   pnew = zeros(size(p));
   pnew(:,1,:) = L*p(:,1,:);
   pnew(:,2,:) = H*p(:,2,:).*(pnew(:,1,:)<=x1 | pnew(:,1,:)>=x2) + ...
                 (p(:,2,:).*(H) + ...
                 (1-p(:,2,:)).*(sqrt(abs(rb^2-(pnew(:,1,:)-L/2).^2)))).* ...
                 (pnew(:,1,:)>x1 & pnew(:,1,:)<x2);
   mesh.dgnodes = pnew;
   mesh.fcurved = repmat(true,size(mesh.f,1),1);
   mesh.tcurved = repmat(true,size(mesh.t,1),1);
end


