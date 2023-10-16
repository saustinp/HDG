function [u,q,uhat,v,rho,iter] = hdgsolve_mah(master, mesh, mesh1, UDG, href, dist, beta, u)

  % href = 1e-3
  % dist = []

if nargin<7
  beta = 1;
end
  
if isempty(mesh1)
  rho = sensor(mesh, master, UDG, href, dist, beta);
else
  rho1 = sensor(mesh1, master, UDG, href, dist, beta);
  [elist, xi] = locatexinmesh(mesh1, mesh.dgnodes(:,1:2,:), [], 1e-4);
  rho = evalfield(mesh1, rho1, elist, xi);
  rho = reshape(rho, [master.nge 1 mesh.ne]); 
end
figure(1); clf; scaplot(mesh, rho);

x = (mesh.dgnodes(:,1,:));
y = (mesh.dgnodes(:,2,:));
L = averagevector(master,mesh);
theta = sum(L(:).*rho(:))/sum(L(:));
if nargin < 8
u = (x.^2+y.^2)/2;
end
uavg = sum(L(:).*u(:));
u = u - uavg;
uhat = inituhat(master,mesh.elcon,u,1);
uhat = uhat(:);

appma.tau = 1;
appma.fbou = 'uboumacyl';
appma.source = 'sourcema5';
appma.neumman = 1;
appma.bcm = [1;1;1;1];
appma.bcd = [1 1; 1 1; 1 1; 1 1];
appma.param = [theta 0 0 0];
[u,q,uhat,v,iter] = hdg_mah(master, mesh, appma, rho, u, uhat);
