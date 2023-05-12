function dist = meshdist(mesh,ib)
% compute the distance to the wall

nd = mesh.nd;
perm = mesh.perm(:,:,1);
ind  = find(mesh.f(:,end)==-ib);

% compute the nodes on the boundary ib
nfe = size(perm,2);
p = [];
for ii=ind'
  k = mesh.f(ii,end-1);  
  fc = mesh.bf(:,k);
  %bf = (fc<0);  
  for j=1:nfe
      if fc(j)==-ib
        p = [p; mesh.dgnodes(perm(:,j),1:nd,k)];        
      end
  end
end

nn = size(mesh.dgnodes,1);
ne = size(mesh.dgnodes,3);
dist = zeros(nn,1,ne);
for i = 1:ne
    for j=1:nn
        x = mesh.dgnodes(j,1:nd,i);
        s = sqrt((p(:,1)-x(1)).^2 + (p(:,2)-x(2)).^2);
        dist(j,1,i) = min(s);  
    end
end
