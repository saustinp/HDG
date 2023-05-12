function udge = enthalpy(pg)

udge = 0*pg(:,1);

gam = 1.4;
Minf = 3.0;
pinf = 1/(gam*Minf^2);
ui = 0.5+pinf/(gam-1);
H = gam*ui - (gam-1)*0.5;
udge(:,1) = H;

