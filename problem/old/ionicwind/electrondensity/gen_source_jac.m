syms ne np nn Er Ez r
syms alpha eta mue mue_ref r_tip k_ep k_np epsilon0 e

f = @(x) x.*(atan(1000*x)/pi + 1/2) + 1/2 - atan(1000)/pi;
d = @(x) atan(1000*x)/pi +(1000*x)/(pi*(1000000*x^2 + 1)) + 1/2;

ne = f(ne);
np = f(np);
nn = f(nn);
% dne_dr = d(dne_dr);
% dne_dz = d(dne_dz);
% dnn_dr = d(dnn_dr);
% dnn_dz = d(dnn_dz);
% dnp_dr = d(dnp_dr);
% dnp_dz = d(dnp_dz);

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

source_electrons = r*((alpha-eta)*(mue/mue_ref)*r_tip*ne*normE - k_ep*epsilon0/(e*mue_ref)*ne*np);
source_negatives = r*(eta*(mue/mue_ref)*r_tip*ne*normE - k_np*epsilon0/(e*mue_ref)*nn*np);
source_positives = r*(alpha*(mue/mue_ref)*r_tip*ne*normE - np*epsilon0/(e*mue_ref)*(k_np*nn + k_ep*ne));
source_poisson = r*(np - ne - nn);

% Note the u vector is: (u1, u2, u3, u4, q1x, q2x, q3x, q4x, q1y, q2y, q3y, q4y)
%                        1   2   3    4   5   6   7   8   9   10  11  12
% Note the u vector is: (ne, np, nn, phi, --, --, --, Er, --, --, --, Ez)
% For the source terms, q[123][xy] are unused hence their derivatives are 0 and are not computed

% Jacobian of electron source terms
dse_dne = diff(source_electrons, 'ne');
dse_dnn = diff(source_electrons, 'nn');
dse_dnp = diff(source_electrons, 'np');
dse_dEr = diff(source_electrons, 'Er');
dse_dEz = diff(source_electrons, 'Ez');

% Jacobian of negative ion source terms
dsn_dne = diff(source_negatives, 'ne');
dsn_dnn = diff(source_negatives, 'nn');
dsn_dnp = diff(source_negatives, 'np');
dsn_dEr = diff(source_negatives, 'Er');
dsn_dEz = diff(source_negatives, 'Ez');

% Jacobian of positive ion source terms
dsp_dne = diff(source_positives, 'ne');
dsp_dnn = diff(source_positives, 'nn');
dsp_dnp = diff(source_positives, 'np');
dsp_dEr = diff(source_positives, 'Er');
dsp_dEz = diff(source_positives, 'Ez');

% Jacobian of electrostatic source terms
dsphi_dne = diff(source_poisson, 'ne');
dsphi_dnn = diff(source_poisson, 'nn');
dsphi_dnp = diff(source_poisson, 'np');

% Writing them all to a file. I then follow it up by a find/replace in the .txt file to replace all the '*' with '.*' etc for the elementwise operations
fid = fopen('source_jac.txt', 'wt');
fprintf(fid, '%s\t', 'dse_dne');
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


