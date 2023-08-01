function [stot, se, sn, sp, sphi] = source_eval(p,udg,param,time)

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

mu_p = param{1};
mu_n = param{2};
Kep = param{5};
Knp = param{6};
e = param{11};
epsilon0 = param{12};
gamma = param{14};
E_bd = param{15};
r_tip = param{16};
n_ref = param{17};
N = param{18};
mue_ref = param{19};
mue = get_mue(normE.*E_bd, N);           % Multiplying by E_bd required to convert back to dimensional units
alpha = get_alpha(normE.*E_bd, N);
eta = get_eta(normE.*E_bd, N);

se = (alpha-eta).*(mue./mue_ref).*r_tip.*ne.*normE - Kep.*epsilon0./(e.*mue_ref).*ne.*np;
sn =         eta.*(mue./mue_ref).*r_tip.*ne.*normE - Knp.*epsilon0./(e.*mue_ref).*nn.*np;
sp =       alpha.*(mue./mue_ref).*r_tip.*ne.*normE -  np.*epsilon0./(e.*mue_ref).*(Knp.*nn + Kep.*ne);
sphi = np - ne - nn;

% Multiply source by r for axisymmetry. I'm leaving the variable name as "f" to avoid changing the code below.
stot = r.*[se sn sp sphi];
