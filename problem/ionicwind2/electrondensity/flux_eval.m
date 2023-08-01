function [f, cr_e, cz_e, cr_n, cz_n, cr_p, cz_p, De_s, Dn_s, Dp_s] = flux_eval(p,udg,param,time)

% Read in values from the p vector -> this normally only contains the (r,z) coordinates but we can also pass in additional fields as appended elements
r = p(:,1,:);
Er0 = p(:,3,:);       % Reading in the values of the E field under the homogeneous forcing (0 space charge) case. Keep in mind that E = -grad(phi), which is what is loaded into p(3,4)
Ez0 = p(:,4,:);

% Read in values from the u vector
ne = udg(:,1,:);
nn = udg(:,2,:);
np = udg(:,3,:);
dne_dr = udg(:,5,:); % q is -grad(u)
dnn_dr = udg(:,6,:);
dnp_dr = udg(:,7,:);
Er_prime = udg(:,8,:);
dne_dz = udg(:,9,:);
dnn_dz = udg(:,10,:);
dnp_dz = udg(:,11,:);
Ez_prime = udg(:,12,:);

Er = Er_prime + Er0;
Ez = Ez_prime + Ez0;
normE = sqrt(Er.^2 + Ez.^2);

% Loading physics parameters
% The electron diffusion and mobility coefficients are functions of the reduced E field
mup = param{1};
mun = param{2};
Dn = param{3};
Dp = param{4};
E_bd = param{15};
r_tip = param{16};
N = param{18};
mue_ref = param{19};
mue = get_mue(normE.*E_bd, N);           % Multiplying by E_bd required to convert back to dimensional units
De = get_diffusion_e(normE.*E_bd, N);

De_s = De./(mue_ref.*E_bd.*r_tip);     % _s is "starred" or nondimensionalized quantity
Dn_s = Dn./(mue_ref.*E_bd.*r_tip);
Dp_s = Dp./(mue_ref.*E_bd.*r_tip);

% Compute convective velocities
cr_e = -(mue./mue_ref).*Er;
cz_e = -(mue./mue_ref).*Ez;
cr_n = -(mun./mue_ref).*Er;
cz_n = -(mun./mue_ref).*Ez;
cr_p = (mup./mue_ref).*Er;
cz_p = (mup./mue_ref).*Ez;

fv = [De_s.*dne_dr, Dn_s.*dnn_dr, Dp_s.*dnp_dr, 0.*Er,...
      De_s.*dne_dz, Dn_s.*dnn_dz, Dp_s.*dnp_dz, 0.*Er];

fi = [cr_e.*ne, cr_n.*nn, cr_p.*np, Er_prime,...
      cz_e.*ne, cz_n.*nn, cz_p.*np, Ez_prime];

% Multiply flux by r for axisymmetry
f = r.*(fi + fv);

end

