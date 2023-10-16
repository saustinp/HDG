function [UDG1, UH1, ACG1, mine1, minf1, ming1, mesh1] = hdgsolve_avloop2(master, mesh, app, q, UDG, UH, lambda0, kappa0, m)
if nargin < 9
    m = 8;
end
mesh1 = mesh;
mesht = mkcgmesh(mesh);
mesht.ib = [];
mesht.in = 1:size(mesht.p2,1);
for d = 1:2
  s = cgpoisson(mesht, master, q(:,d,:), [0 1.0]);    
  s = s(mesht.t2');
  mesh1.dgnodes(:,d,:) = reshape(s, size(q(:,d,:)));
end
mesh1.dgnodes(:,3,:) = 0;

S0 = 0.2;
eta = 0.9; %m = 8;
lambda = ones(m,1)*lambda0;
for i = 2:m
    lambda(i) = lambda(i-1)*eta;
end
kappa = ones(m,1)*kappa0;
for i = 2:m
    kappa(i) = 1 + (kappa(i-1)-1)*eta;
end

[UDG1, UH1, ACG1, mine1, minf1, ming1] = avloop(master, mesh1, app, UDG, UH, S0, lambda, kappa);
