% These inputs must match the dimensions of the input vectors in the final output function. For example:
% function [s,jac_out] = flux_electrondensity(p,udg,param,time)

syms x1 x2 x3 x4  % Variables stored in the p vector
syms u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12     % Components of UDG
syms time
syms param1 param2 param3 param4 param5 param6 param7 param8 param9 param10 param11 param12 param13 param14 param15 param16 param17 param18 param19     % Components of the physics parameters

param = [param1 param2 param3 param4 param5 param6 param7 param8 param9 param10 param11 param12 param13 param14 param15 param16 param17 param18 param19];
udg = [u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12];
p = [x1 x2 x3 x4];

nd = 2;
ncu = 4;    % num components of U
nch = ncu;    % num components of UHAT
nc = 12;    % num components of UDG

% Read in values from the p vector -> this normally only contains the (r,z) coordinates but we can also pass in additional fields as appended elements
r = p(1);
Er0 = p(3);       % Reading in the values of the E field under the homogeneous forcing (0 space charge) case. Keep in mind that E = -grad(phi), which is what is loaded into p(3,4)
Ez0 = p(4);

% Read in values from the u vector
ne = udg(1);
nn = udg(2);
np = udg(3);

% clip_fcn = @(x) x*(atan(1000*x)/pi + 1/2) + 1/2 - atan(1000)/pi;
% ne = clip_fcn(udg(1));
% nn = clip_fcn(udg(2));
% np = clip_fcn(udg(3));

dne_dr = udg(5); % q is -grad(u)
dnn_dr = udg(6);
dnp_dr = udg(7);
Er_prime = udg(8);
dne_dz = udg(9);
dnn_dz = udg(10);
dnp_dz = udg(11);
Ez_prime = udg(12);

Er = Er_prime + Er0;
Ez = Ez_prime + Ez0;
normE = sqrt(Er.^2 + Ez.^2);

% Loading physics parameters
% The electron diffusion and mobility coefficients are functions of the reduced E field
mup = param(1);
mun = param(2);
Dn = param(3);
Dp = param(4);
Kep = param(5);
Knp = param(6);
Nmax = param(7);
r0 = param(8);
z0 = param(9);
s0 = param(10);
e = param(11);
epsilon0 = param(12);
Ua = param(13);
gamma = param(14);
E_bd = param(15);
r_tip = param(16);
n_ref = param(17);
N = param(18);
mue_ref = param(19);

De = get_diffusion_e(normE*E_bd, N)*10;
mue = get_mue(normE*E_bd, N);           % Multiplying by E_bd required to convert back to dimensional units
alpha = get_alpha(normE*E_bd, N);
eta = get_eta(normE*E_bd, N);


%%%% FLUX %%%%
% De_s = De/(mue_ref*E_bd*r_tip);     % _s is "starred" or nondimensionalized quantity
% Dn_s = Dn/(mue_ref*E_bd*r_tip);
% Dp_s = Dp/(mue_ref*E_bd*r_tip);

% Compute convective velocities
cr_e = -mue.*Er;
cz_e = -mue.*Ez;
cr_n = -mun.*Er;
cz_n = -mun.*Ez;
cr_p = mup.*Er;
cz_p = mup.*Ez;

fv = [De.*dne_dr, Dn.*dnn_dr, Dp.*dnp_dr, -Er_prime,...
      De.*dne_dz, Dn.*dnn_dz, Dp.*dnp_dz, -Ez_prime];

fi = [cr_e.*ne, cr_n.*nn, cr_p.*np, 0,...
      cz_e.*ne, cz_n.*nn, cz_p.*np, 0];

% Multiply flux by r for axisymmetry
f = r*(fi + fv);


%%%% SOURCE %%%%
se = (alpha-eta)*mue*ne*normE - Kep*ne*np;
sn =         eta*mue*ne*normE - Knp*nn*np;
sp =       alpha*mue*ne*normE - np*(Knp*nn + Kep*ne);
sphi = ne + nn - np;          % Because we have to enter the -grad(u) for the flux

% Multiply source by r for axisymmetry. I'm leaving the variable name as "s" to avoid changing the code below.
s = r*[se sn sp sphi];

gencode(p, f, udg, param, time, 'flux', 2, 4);
disp('Wrote flux... done')
gencode(p, s, udg, param, time, 'source', 2, 4);
disp('Wrote source... done')

check_gencode;