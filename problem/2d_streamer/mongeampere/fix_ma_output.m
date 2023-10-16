% mesh1 = mesh2;      % mesh2 is the undeformed mesh pre-MA
% mesht = mkcgmesh(mesh2);
% mesht.ib = [];
% mesht.in = 1:size(mesht.p2,1);
% for d = 1:2
%   s = cgpoisson(mesht, master, q(:,d,:), [0 1.0]);    
%   s = s(mesht.t2');
%   mesh1.dgnodes(:,d,:) = reshape(s, size(q(:,d,:)));
% end

[cgnodes,cgelcon,rowent2elem,colent2elem,cgent2dgent] = mkcgent2dgent(mesh2.dgnodes(:,1:2,:),1e-8);
q = ucg2udg(udg2ucg(q, cgent2dgent, rowent2elem), cgelcon);
mesh1.dgnodes = q;
meshplot(mesh1);