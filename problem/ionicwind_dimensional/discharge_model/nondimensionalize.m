syms ne nn np phi E t r z
syms ne_s nn_s np_s phi_s E_s t_s r_s z_s
syms mue mun mup kep knp De Dn Dp alpha eta
syms epsilon0 E_bd e r_tip mue_ref

ne = ne_s*epsilon0*E_bd/(e*r_tip);
nn = nn_s*epsilon0*E_bd/(e*r_tip);
np = np_s*epsilon0*E_bd/(e*r_tip);
E = E_s*E_bd;
phi = phi_s*E_bd*r_tip;
t = t_s*r_tip/(mue_ref*E_bd);

factor_drift = ne/t / (ne_s/t_s);
factor_poisson = phi/r_tip^2;

% Don't forget to include a 1/l_ref for each gradient term
% eq1t1 = ne/t /factor_drift
% eq1t21 = (- De*ne/r_tip)/r_tip /factor_drift
% eq1t22 = (-mue*E*ne)/r_tip /factor_drift
% eq1t3 = (alpha-eta)*ne*mue*abs(E) /factor_drift
% eq1t4 = - kep*ne*np /factor_drift

% eq2t1 = nn/t /factor_drift
% eq2t21 = (-mun*E*nn)/r_tip /factor_drift
% eq2t22 = (- Dn*nn/r_tip)/r_tip /factor_drift
% eq2t3 = eta*ne*mue*abs(E) /factor_drift
% eq2t4 = - knp*nn*np /factor_drift

% eq3t1 = np/t /factor_drift
% eq3t21 = (mup*E*np)/r_tip /factor_drift
% eq3t22 = - (Dp*np/r_tip)/r_tip /factor_drift
% eq3t3 = alpha*ne*mue*abs(E) /factor_drift
% eq3t4 = - knp*nn*np /factor_drift
% eq3t5 = - kep*ne*np /factor_drift

eqn4 = phi/r_tip^2 == -e/epsilon0*(np-ne-nn);
eqn4_nondim = solve(eqn4, phi_s)

% eq1 = eq1t1+eq1t2==eq1t3+eq1t4
% eq2 = eq2t1+eq2t2==eq2t3+eq2t4
% eq3 = eq3t1+eq3t2==eq3t3+eq3t4+eq3t5


% Output
% eq1t1 = ne_s/t_s
% eq1t21 = -(De*ne_s)/(E_bd*mue_ref*r_tip) 
% eq1t22 = -(E_s*mue*ne_s)/mue_ref
% eq1t3 = (mue*ne_s*r_tip*abs(E_bd*E_s)*(alpha - eta))/(E_bd*mue_ref)
% eq1t4 = -(epsilon0*kep*ne_s*np_s)/(e*mue_ref)

% eq2t1 = nn_s/t_s
% eq2t2 = -(e*r_tip*((Dn*E_bd*epsilon0*nn_s)/(e*r_tip^2) + (E_bd^2*E_s*epsilon0*mun*nn_s)/(e*r_tip)))/(E_bd^2*epsilon0*mue_ref)
% eq2t3 = (eta*mue*ne_s*r_tip*abs(E_bd*E_s))/(E_bd*mue_ref)
% eq2t4 = -(epsilon0*knp*nn_s*np_s)/(e*mue_ref)

% eq3t1 = np_s/t_s
% eq3t21 = (E_s*mup*np_s)/mue_ref
% eq3t22 = -(Dp*np_s)/(E_bd*mue_ref*r_tip)
% eq3t3 = (alpha*mue*ne_s*r_tip*abs(E_bd*E_s))/(E_bd*mue_ref)
% eq3t4 = -(epsilon0*knp*nn_s*np_s)/(e*mue_ref)
% eq3t5 = -(epsilon0*kep*ne_s*np_s)/(e*mue_ref)