
% Mach   = 8.03;
% aoa    = 0.0;
% gam = 1.4;
% epslm = 0.0;
% Minf = Mach;                  % Infinity conditions
% pinf = 1/(gam*Minf^2);
% Re = 1.94e5;
% Pr = 0.71;
% alpha = aoa*pi/180;
% tau = 4;
% kk = 30;
% Tref  = 122.11;
% Twall = 294.44;
param = init_phys_param();




n = 1;
pg = rand(n,4);
udg = rand(n,12);
[f,f_udg] = source2d(pg,udg,param,0);

g_udg = 0*f_udg;
epsilon = 1e-6;
for i = 1:12
  u1 = udg;
  u1(i) = u1(i) + epsilon;
  [f1,t1] = source2d(pg,u1,param,0);
  
  u2 = udg;
  u2(i) = u2(i) - epsilon;
  [f2,t2] = source2d(pg,u2,param,0);
  g_udg(:,:,i) = (f1-f2)/(2*epsilon);
end
max(abs(f_udg(:) - g_udg(:)))

