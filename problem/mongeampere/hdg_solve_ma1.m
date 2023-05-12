function [udg,uhat] = hdg_solve_ma1(master, mesh, app, udg)
% Extend hdg_ma2 to deal with neumman boundary condition

[Ae, ~, De, ~, Ue, MiB, MiC, MiD, MiE, G, L, H, ~, P] = hdg_precompute(master, mesh, app, udg);

nd = master.nd;
vold = udg(:,2:(nd+1),:);
qold = udg(:,(nd+2):end,:);

itermax = 60;
iter = 0;
while (iter<itermax)
  iter = iter + 1;  
  
  F = hdg_sourcevector(master, app, dgnodes, udg);
  [S, G] = hdg_boundaryvector(master, app, mesh.dgnodes, udg, G);
  [Ae, Re, We] = hdg_assemble(master, mesh, app, Ae, MiB, MiC, De, Ue, F, H, G, L, S);

  [c, uc] = uaverage(We, Ue, P, mesh.elcon);  

  uhat = hdg_linearsystem(Ae, Re, mesh.elcon, app.neumman, uc, c);  

  u = getu(We, Ue, uhat, mesh.elcon);
  q = getq(MiB, MiC, udg(:,1,:), uhat, mesh.elcon);
  v = getv(MiD, MiE, udg(:,1,:), uhat, mesh.elcon);

  udg(:,1,:) = u;
  udg(:,2:(nd+1),:) = q;
  udg(:,(nd+2):(1+nd+nd*nd),:) = v;
    
  fprintf('Iteration: %d,  ||q - qold|| = %e,   ||v - vold|| = %e\n', [iter norm(q(:)-qold(:)) norm(v(:)-vold(:))]);         
  if norm(v(:)-vold(:)) <= 1e-4
    break;
  end
  
  vold = v;
  qold = q;  
end


end

