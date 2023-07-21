% Script to check the jacobian of the flux and source functions using a finite difference
param = init_phys_param();

nc_u = 4;
ndim = 2;
nc_udg = nc_u*(ndim+1);
n_eval_pts = 1;       % Number of points to evalauate the flux and source functions at

pg = rand(n_eval_pts,nc_u);
udg = rand(n_eval_pts,nc_udg);
[~,f_udg] = flux2d(pg,udg,param,0);
[~,s_udg] = source2d(pg,udg,param,0);

f_udg_FD = 0*f_udg;   % Initialize empty arrays to hold the finite differenced source and flux jacobians
s_udg_FD = 0*s_udg;

eps = 1e-6;
for i = 1:nc_udg
  u1 = udg;
  u2 = udg;
  u1(i) = u1(i) + eps;      % Perturb one of the components of u by +eps
  u2(i) = u2(i) - eps;      % Perturb the same component of u by -eps

  [f1,~] = flux2d(pg,u1,param,0);
  [s1,~] = source2d(pg,u1,param,0);
  
  [f2,~] = flux2d(pg,u2,param,0);
  [s2,~] = source2d(pg,u2,param,0);

  f_udg_FD(:,:,:,i) = (f1-f2)/(2*eps);    % 2nd order FD
  s_udg_FD(:,:,i) = (s1-s2)/(2*eps);
end

flux_max_error = max(abs(f_udg(:) - f_udg_FD(:)))
source_max_error = max(abs(s_udg(:) - s_udg_FD(:)))