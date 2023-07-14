% These inputs must match the dimensions of the input vectors in the final output function. For example:
% function [f,f_udg] = flux_electrondensity(p,udg,param,time)

syms x1 x2 x3 x4  % Variables stored in the p vector
syms u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12     % Components of UDG
syms time
syms param1 param2 param3 param4 param5 param6 param7 param8 param9 param10 param11 param12 param13 param14 param15 param16 param17 param18 param19     % Components of the physics parameters
syms zero one

param = [param1 param2 param3 param4 param5 param6 param7 param8 param9 param10 param11 param12 param13 param14 param15 param16 param17 param18 param19];
udg = [u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12];
p = [x1 x2 x3 x4];

% Read in values from the p vector -> this normally only contains the (r,z) coordinates but we can also pass in additional fields as appended elements
r = p(1);
Er0 = p(3);       % Reading in the values of the E field under the homogeneous forcing (0 space charge) case. Keep in mind that E = -grad(phi), which is what is loaded into p(3,4)
Ez0 = p(4);

% Read in values from the u vector
ne = udg(1);
nn = udg(2);
np = udg(3);
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
mu_p = param(1);
mu_n = param(2);
Kep = param(5);
Knp = param(6);
e = param(11);
epsilon0 = param(12);
gamma = param(14);
E_bd = param(15);
r_tip = param(16);
n_ref = param(17);
N = param(18);
mue_ref = param(19);
mue = get_mue(normE*E_bd, N);           % Multiplying by E_bd required to convert back to dimensional units
alpha = get_alpha(normE*E_bd, N);
eta = get_eta(normE*E_bd, N);

se = (alpha-eta)*(mue/mue_ref)*r_tip*ne*normE - Kep*epsilon0/(e*mue_ref)*ne*np;
sn =         eta*(mue/mue_ref)*r_tip*ne*normE - Knp*epsilon0/(e*mue_ref)*nn*np;
sp =       alpha*(mue/mue_ref)*r_tip*ne*normE -  np*epsilon0/(e*mue_ref)*(Knp*nn + Kep*ne);
sphi = np - ne - nn;

% Multiply source by r for axisymmetry. I'm leaving the variable name as "f" to avoid changing the code below.
f1 = r*[se sn sp sphi];

nd = 2;
ncu = 4;    % num components of U
nch = ncu;    % num components of UHAT
nc = 12;    % num components of UDG

for n=1:1
    if n==1
        f = f1;
        filename1 = ['source' num2str(nd) 'd' '.m'];
    elseif n==2
        f = f2;
        filename1 = ['source' num2str(nd) 'd2' '.m'];
    end
    
    %%% compute Jacobian
    jac_f  = jacobian(f,udg);

    %%% And patch with vector zero or one to have the right sizes
    for ii = 1:size(f,1)
        for jj = 1:size(f,2)
            temp = ismember(symvar(f(ii,jj)),udg);
            if f(ii,jj)==0, f(ii,jj) = zero; %%% No dependency on anything
            elseif isempty(temp) || sum(temp)==0, f(ii,jj) = (f(ii,jj))*one; %%% No dependency on state vars
            end
        end
    end

    for ii = 1:size(jac_f,1)
        for jj = 1:size(jac_f,2)
            temp = ismember(symvar(jac_f(ii,jj)),udg);
            if jac_f(ii,jj)==0, jac_f(ii,jj) = zero; %%% No dependency on anything
            elseif isempty(temp) || sum(temp)==0, jac_f(ii,jj) = (jac_f(ii,jj))*one; %%% No dependency on state vars
            end
        end
    end    
    
    %matlabFunction(f(1),'file','tmp.m','vars',{pg,udg,param,time,[zero one]},'outputs', {'f'});
    
    % generate a temporary matlab file         $ NOTE: changed "pg" -> "p"
    matlabFunction(f(:),jac_f(:),'file','tmp.m','vars',{p,udg,param,time,[zero one]},'outputs', {'s','s_udg'});          % NOTE: This was changed from the template

    % open the file and modify it
    fid = fopen('tmp.m','r');
    gid = fopen(filename1,'wt');

    tline = fgetl(fid); 
    i=1;       
    while ischar(tline)        
        str = strrep(tline, 'tmp', strrep(filename1,'.m',''));
        str = strrep(str, 'TMP', upper(strrep(filename1,'.m','')));
        str = strrep(str, 'in1', 'p');          % NOTE: This was changed from the template
        str = strrep(str, 'IN1', 'P');        
        str = strrep(str, 'in2', 'udg');
        str = strrep(str, 'IN2', 'UDG');        
        str = strrep(str, 'in3', 'param');
        str = strrep(str, 'IN3', 'PARAM');                
        str = strrep(str, ',in5)', ')');                
        str = strrep(str, ',IN5)', ')');    
        str = strrep(str, 'param(:,1)', 'param{1}'); 
        str = strrep(str, 'param(:,2)', 'param{2}'); 
        str = strrep(str, 'param(:,3)', 'param{3}'); 
        str = strrep(str, 'param(:,4)', 'param{4}'); 
        str = strrep(str, 'param(:,5)', 'param{5}'); 
        str = strrep(str, 'param(:,6)', 'param{6}'); 
        str = strrep(str, 'param(:,7)', 'param{7}'); 
        str = strrep(str, 'param(:,8)', 'param{8}'); 
        str = strrep(str, 'param(:,9)', 'param{9}');     
        str = strrep(str, 'param(:,10)', 'param{10}');     
        str = strrep(str, 'param(:,11)', 'param{11}');       % NOTE: These lines were added
        str = strrep(str, 'param(:,12)', 'param{12}');     
        str = strrep(str, 'param(:,13)', 'param{13}');     
        str = strrep(str, 'param(:,14)', 'param{14}');     
        str = strrep(str, 'param(:,15)', 'param{15}');     
        str = strrep(str, 'param(:,16)', 'param{16}');     
        str = strrep(str, 'param(:,17)', 'param{17}');     
        str = strrep(str, 'param(:,18)', 'param{18}');     
        str = strrep(str, 'param(:,19)', 'param{19}');     
        str = strrep(str, 'in5(:,1)', 'zeros(ng,1)');
        str = strrep(str, 'in5(:,2)', 'ones(ng,1)');        
        if i==7
            str = '[ng,nc] = size(udg);';
            fprintf(gid, '%s\n', str);                  
            str = ['nch = ' num2str(nch) ';'];
            fprintf(gid, '%s\n', str);                  
            str = ['nd = ' num2str(nd) ';'];
        end
        fprintf(gid, '%s\n', str);                  
        tline = fgetl(fid);        
        i=i+1;
        %disp(str)
    end

    % NOTE: This will produce a function that has errors. You need to move the following "reshape" lies above the second "end" statement in order to be enclosed in the the function
    str = 's = reshape(s,ng,nch);';          % NOTE: This was changed from the template
    fprintf(gid, '%s\n', str);                  
    str = 's_udg = reshape(s_udg,ng,nch,nc);';          % NOTE: This was changed from the template
    fprintf(gid, '%s\n', str);                  

    fclose(fid);
    fclose(gid);
    delete('tmp.m');
end
