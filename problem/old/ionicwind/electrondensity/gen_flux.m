syms ne np nn Er Ez r dne_dr dne_dz dnn_dr dnn_dz dnp_dr dnp_dz
syms alpha eta mue mue_ref r_tip k_ep k_np epsilon0 e
syms cr_e cr_n cr_p cz_e cz_n cz_p Ex_prime Ey_prime De_star Dn_star Dp_star

% f = @(x) x.*(atan(1000*x)/pi + 1/2) + 1/2 - atan(1000)/pi;
% d = @(x) atan(1000*x)/pi +(1000*x)/(pi*(1000000*x^2 + 1)) + 1/2;

% ne = f(ne);
% np = f(np);
% nn = f(nn);
% % dne_dr = d(dne_dr);
% % dne_dz = d(dne_dz);
% % dnn_dr = d(dnn_dr);
% % dnn_dz = d(dnn_dz);
% % dnp_dr = d(dnp_dr);
% % dnp_dz = d(dnp_dz);

N = 2.4614924955148245e+25;     % Neutral number density
normE = sqrt(Er^2+Ez^2);
swarm = get_swarm_params(normE, N);
alpha = swarm(:,1);
eta = swarm(:,2);
beta = swarm(:,3);
De = swarm(:,4);
mue = swarm(:,5);
mup = swarm(:,6);
mun = swarm(:,7);

% source_electrons = r*((alpha-eta)*(mue/mue_ref)*r_tip*ne*normE - k_ep*epsilon0/(e*mue_ref)*ne*np);
% source_negatives = r*(eta*(mue/mue_ref)*r_tip*ne*normE - k_np*epsilon0/(e*mue_ref)*nn*np);
% source_positives = r*(alpha*(mue/mue_ref)*r_tip*ne*normE - np*epsilon0/(e*mue_ref)*(k_np*nn + k_ep*ne));  
% source_poisson = r*(np - ne - nn);

% Fluxes
fe_r = cr_e.*(r.*ne) + De_star.*(r.*dne_dr);
fe_z = cz_e.*(r.*ne) + De_star.*(r.*dne_dz);
fn_r = cr_n.*(r.*nn) + Dn_star.*(r.*dnn_dr);
fn_z = cz_n.*(r.*nn) + Dn_star.*(r.*dnn_dz);
fp_r = cr_p.*(r.*np) + Dp_star.*(r.*dnp_dr);
fp_z = cz_p.*(r.*np) + Dp_star.*(r.*dnp_dz);
fphi_r = r.*Ex_prime;
fphi_z = r.*Ey_prime;

% Note the u vector is: (u1, u2, u3, u4, q1x, q2x, q3x, q4x, q1y, q2y, q3y, q4y)
%                        1   2   3    4   5   6   7   8   9   10  11  12
% Note the u vector is: (ne, np, nn, phi, --, --, --, Er, --, --, --, Ez)
% For the source terms, q[123][xy] are unused hence their derivatives are 0 and are not computed



v11 = diff(fe_r, 'ne')
v12 = diff(fe_r, 'dne_dr')
v13 = diff(fe_z, 'ne')
v14 = diff(fe_z, 'dne_dz')
v15 = diff(fe_r, 'Er')
v16 = diff(fe_z, 'Ez')

v21 = diff(fn_r, 'nn')
v22 = diff(fn_r, 'dnn_dr')
v23 = diff(fn_z, 'nn')
v24 = diff(fn_z, 'dnn_dz')
v25 = diff(fn_r, 'Er')
v22 = diff(fn_z, 'Ez')

v31 = diff(fp_r, 'np')
v32 = diff(fp_r, 'dnp_dr')
v33 = diff(fp_z, 'np')
v34 = diff(fp_z, 'dnp_dz')
v35 = diff(fp_r, 'Er')
v36 = diff(fp_z, 'Ez')

% v41 = diff(fphi_r, 'Er')
% v42 = diff()







% % Jacobian of electron source terms
% dse_dne = diff(source_electrons, 'ne');
% dse_dnn = diff(source_electrons, 'nn');
% dse_dnp = diff(source_electrons, 'np');
% dse_dEr = diff(source_electrons, 'Er');
% dse_dEz = diff(source_electrons, 'Ez');

% % Jacobian of negative ion source terms
% dsn_dne = diff(source_negatives, 'ne');
% dsn_dnn = diff(source_negatives, 'nn');
% dsn_dnp = diff(source_negatives, 'np');
% dsn_dEr = diff(source_negatives, 'Er');
% dsn_dEz = diff(source_negatives, 'Ez');

% % Jacobian of positive ion source terms
% dsp_dne = diff(source_positives, 'ne');
% dsp_dnn = diff(source_positives, 'nn');
% dsp_dnp = diff(source_positives, 'np');
% dsp_dEr = diff(source_positives, 'Er');
% dsp_dEz = diff(source_positives, 'Ez');

% % Jacobian of electrostatic source terms
% dsphi_dne = diff(source_poisson, 'ne');
% dsphi_dnn = diff(source_poisson, 'nn');
% dsphi_dnp = diff(source_poisson, 'np');

% Writing them all to a file. I then follow it up by a find/replace in the .txt file to replace all the '*' with '.*' etc for the elementwise operations
fid = fopen('source_jac.txt', 'wt');
fprintf(fid, '%s\t', 'v11');
fprintf(fid, '%s\n', char(dse_dne));
fprintf(fid, '%s\t', 'dse_dnn');
fprintf(fid, '%s\n', char(dse_dnn));
fprintf(fid, '%s\t', 'dse_dnp');
fprintf(fid, '%s\n', char(dse_dnp));
fprintf(fid, '%s\t', 'dse_dEr');
fprintf(fid, '%s\n', char(dse_dEr));
fprintf(fid, '%s\t', 'dse_dEz');
fprintf(fid, '%s\n', char(dse_dEz));

fprintf(fid, '%s\t', 'dsn_dne');
fprintf(fid, '%s\n', char(dsn_dne));
fprintf(fid, '%s\t', 'dsn_dnn');
fprintf(fid, '%s\n', char(dsn_dnn));
fprintf(fid, '%s\t', 'dsn_dnp');
fprintf(fid, '%s\n', char(dsn_dnp));
fprintf(fid, '%s\t', 'dsn_dEr');
fprintf(fid, '%s\n', char(dsn_dEr));
fprintf(fid, '%s\t', 'dsn_dEz');
fprintf(fid, '%s\n', char(dsn_dEz));

fprintf(fid, '%s\t', 'dsp_dne');
fprintf(fid, '%s\n', char(dsp_dne));
fprintf(fid, '%s\t', 'dsp_dnn');
fprintf(fid, '%s\n', char(dsp_dnn));
fprintf(fid, '%s\t', 'dsp_dnp');
fprintf(fid, '%s\n', char(dsp_dnp));
fprintf(fid, '%s\t', 'dsp_dEr');
fprintf(fid, '%s\n', char(dsp_dEr));
fprintf(fid, '%s\t', 'dsp_dEz');
fprintf(fid, '%s\n', char(dsp_dEz));

fprintf(fid, '%s\t', 'dsphi_dne');
fprintf(fid, '%s\n', char(dsphi_dne));
fprintf(fid, '%s\t', 'dsphi_dnn');
fprintf(fid, '%s\n', char(dsphi_dnn));
fprintf(fid, '%s\t', 'dsphi_dnp');
fprintf(fid, '%s\n', char(dsphi_dnp));
fclose(fid);


