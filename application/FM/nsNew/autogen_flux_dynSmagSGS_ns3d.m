
% DESCRIPTION:

% Momentum equation. Only extra-term is SGS viscous stress:
% - Deviatoric SGS viscous stress tensor model: Dynamic Smagorinksy
% - Isotropic SGS viscous stress tensor model: Yoshizawa

% Energy equation. Three extre-terms: SGS heat transfer, SGS turbulent
% diffusion, SGS viscous diffusion. The 2 first terms are of similar
% magnitude. The third one is about 1 order of magnitude smaller
% - SGS heat transfer: Turbulent Prandtl number based approach
% - SGS turbulent diffusion: Knight et al. (1998) approach (u_i *
% tau_{SGS,j,i}. This could be the only attempt to model this term.
% Alternatively, it can be modeled together with the SGS heat transfer.
% - SGS viscous diffusion: Neglected

fileName = 'flux_dynSmagSGS_ns3d';

% Input variables:
syms r ru rv rw rE rx rux rvx rwx rEx ry ruy rvy rwy rEy rz ruz rvz rwz rEz
syms ux uy uz vx vy vz wx wy wz T p
syms gam epslm Re Minf Pr tau gam1
syms pg time
syms S_xx S_xy S_xz S_yx S_yy S_yz S_zx S_zy S_zz S_mag nu_SGS
syms h C_d C_i Pr_SGS Delta
syms x y z h

udg = [r ru rv rw rE rx rux rvx rwx rEx ry ruy rvy rwy rEy rz ruz rvz rwz rEz];
pg = [x y z h];
param = [gam epslm Re Pr Minf tau];

% C_i = 0.0066;
% Values for C_i:
% 1) 0.09
% 2) E. Garnier, N. Adams, P. Sagaut, Large Eddy Simulation for Compressible
% Flows : 0.0066-0.066 (recommend 0.0066, although no big differences
% have been observed within this range)
% 2) T.B. Gatski, J.P. Bonnet, Compressibility, Turbulence and High Speed
% Flow: 0.09 (suggested by Yoshizawa), 0.0066 (from calibration with DNS)

% Pr_SGS = 0.5;
% Values for Pr_SGS:
% 1) T.B. Gatski, J.P. Bonnet, Compressibility, Turbulence and High Speed
% Flow: 0.5 (some group assumed this value to calibrate C_d, C_i with DNS data
% 2) M.P. Martin, U. Piomelli, Subgrid-Scale Models for Compressible
% Large-Eddy Simulations: Used a value of 0.7 for numerical experiments and
% comparison with DNS
% 3) CHARLES: 0.9

% Delta = V^(1/3);

% Field variables:
u = ru./r;
v = rv./r;
w = rw./r;
q = 0.5*(u.*u+v.*v+w.*w);
p = gam1*(rE-r.*q);
T = p./(gam1*r);

ux = (rux - u.*rx)./r;
vx = (rvx - v.*rx)./r;
wx = (rwx - w.*rx)./r;
uy = (ruy - u.*ry)./r;
vy = (rvy - v.*ry)./r;
wy = (rwy - w.*ry)./r;
uz = (ruz - u.*rz)./r;
vz = (rvz - v.*rz)./r;
wz = (rwz - w.*rz)./r;

qx = u.*ux + v.*vx + w.*wx;
qy = u.*uy + v.*vy + w.*wy;
qz = u.*uz + v.*vz + w.*wz;

px = gam1*(rEx - rx.*q - r.*qx);
py = gam1*(rEy - ry.*q - r.*qy);
pz = gam1*(rEz - rz.*q - r.*qz);

Tx  = (px.*r - p.*rx)./(gam1*r.^2);
Ty  = (py.*r - p.*ry)./(gam1*r.^2);
Tz  = (pz.*r - p.*rz)./(gam1*r.^2);

S_xx = ux;
S_xy = 0.5 * (vx + uy);
S_xz = 0.5 * (wx + uz);
S_yx = 0.5 * (uy + vx);
S_yy = vy;
S_yz = 0.5 * (wy + vz);
S_zx = 0.5 * (uz + wx);
S_zy = 0.5 * (vz + wy);
S_zz = wz;

S_mag = sqrt( 2.0 * (S_xx*S_xx + S_xy*S_xy + S_xz*S_xz + ...
                     S_yx*S_yx + S_yy*S_yy + S_yz*S_yz + ...
                     S_zx*S_zx + S_zy*S_zy + S_zz*S_zz) );

nu_SGS = (C_d * Delta)^2 * S_mag;

txx_SGS = r.*nu_SGS.*(2*ux - (2./3)*(ux+vy+wz)) + (2/3)*r*C_i*Delta^2*S_mag^2;
txy_SGS = r.*nu_SGS.*(uy+vx);
txz_SGS = r.*nu_SGS.*(uz+wx);
tyx_SGS = r.*nu_SGS.*(vx+uy);
tyy_SGS = r.*nu_SGS.*(2*vy - (2./3)*(ux+vy+wz)) + (2/3)*r*C_i*Delta^2*S_mag^2;
tyz_SGS = r.*nu_SGS.*(vz+wy);
tzx_SGS = r.*nu_SGS.*(wx+uz);
tzy_SGS = r.*nu_SGS.*(wy+vz);
tzz_SGS = r.*nu_SGS.*(2*wz - (2./3)*(ux+vy+wz)) + (2/3)*r*C_i*Delta^2*S_mag^2;

factor1 = (gam*r.*nu_SGS)./Pr_SGS;

F11 = 0;
F12 = 0;
F13 = 0;
F21 = txx_SGS;
F22 = txy_SGS;
F23 = txz_SGS;
F31 = tyx_SGS;
F32 = tyy_SGS;
F33 = tyz_SGS;
F41 = tzx_SGS;
F42 = tzy_SGS;
F43 = tzz_SGS;
F51 = 0.5 * (txx_SGS.*u + txy_SGS.*v + txz_SGS.*w) + factor1.*Tx;
F52 = 0.5 * (tyx_SGS.*u + tyy_SGS.*v + tyz_SGS.*w) + factor1.*Ty;
F53 = 0.5 * (tzx_SGS.*u + tzy_SGS.*v + tzz_SGS.*w) + factor1.*Tz;

F_SGS = [F11; F21; F31; F41; F51; F12; F22; F32; F42; F52; F13; F23; F33; F43; F53];

F_SGS_udg = jacobian(F_SGS,udg);

F_SGS = reshape(F_SGS,[5, 3]);
F_SGS_udg = reshape(F_SGS_udg,[5, 3, 20]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate C code:
disp('  ')
disp('  ')
disp('Auto-generating C code...'); tic

temporaryFileName1 = ['tmp_', fileName, '1.c'];
temporaryFileName2 = ['tmp_', fileName, '2.c'];

% generate temporary C file
ccode(F_SGS(:),'file',temporaryFileName1)
ccode(F_SGS_udg(:),'file',temporaryFileName2)

fileNameC = [fileName, '.c'];
delete(fileNameC);
gid = fopen(fileNameC,'wt');

fprintf(gid, '\n');
str = '// Written by: C. Nguyen & P. Fernandez';
fprintf(gid, '%s\n\n',str);
str = ['void ' fileName '(double *f, double *f_udg, double *pg, double *udg, double *Cs, appstruct &app, double *param, double time, int ng, int nc, int ncu, int nd, int ncd, int computeJacobian)'];
fprintf(gid, '%s\n', str); 
str = '{';
fprintf(gid, '%s\n', str);         

str = 'double gam = param[0];';
fprintf(gid, '\n\t%s\n', str);   
str = 'double gam1 = gam - 1.0;';
fprintf(gid, '\t%s\n', str);    
str = 'double C_i = 0.0066;';
fprintf(gid, '\t%s\n', str);   
str = 'double Pr_SGS = 0.5;';
fprintf(gid, '\t%s\n', str);   

% for i=1:length(param)
%     str = ['double param' num2str(i) ' = param[' num2str(i-1) '];'];
%     fprintf(gid, '\t%s\n', str);                  
% end        

%     str = 'double ';
%     for i=1:length(pg)-1
%         str = [str 'x' num2str(i) ','];        
%     end
%     str = [str 'x' num2str(length(pg)) ';'];        
%     fprintf(gid, '\t%s\n', str);         
%     
%     str = 'double ';
%     for i=1:length(udg)-1
%         str = [str 'u' num2str(i) ','];        
%     end
%     str = [str 'u' num2str(length(udg)) ';'];        
%     fprintf(gid, '\t%s\n', str);         

str = 'for (int i = 0; i <ng; i++) {';
fprintf(gid, '\n\t%s\n', str);         

% for i=1:length(pg)
%     str = ['double x' num2str(i) ' = pg[' num2str(i-1) '*ng+i];'];
%     fprintf(gid, '\t\t%s\n', str);                  
% end
str = 'double x = pg[0*ng+i];';
fprintf(gid, '\t\t%s\n', str);
str = 'double y = pg[1*ng+i];';
fprintf(gid, '\t\t%s\n', str);
str = 'double z = pg[2*ng+i];';
fprintf(gid, '\t\t%s\n', str);
str = 'double h = pg[3*ng+i];';
fprintf(gid, '\t\t%s\n', str);
str = 'double Delta = h;        // TODO: Implement Delta = V^(1/3)';
fprintf(gid, '\t\t%s\n', str);
str = 'double C_d = Cs[i];';
fprintf(gid, '\t\t%s\n', str);
fprintf(gid, '\t\t\n');
for i=1:length(udg)
    str = ['double u' num2str(i) ' = udg[' num2str(i-1) '*ng+i];'];
    fprintf(gid, '\t\t%s\n', str);                  
end
fprintf(gid, '\n');      

str = 'double r = u1;';
fprintf(gid, '\t\t%s\n', str);
str = 'double ru = u2;';
fprintf(gid, '\t\t%s\n', str);
str = 'double rv = u3;';
fprintf(gid, '\t\t%s\n', str);
str = 'double rw = u4;';
fprintf(gid, '\t\t%s\n', str);
str = 'double rE = u5;';
fprintf(gid, '\t\t%s\n', str);
str = 'double rx = u6;';
fprintf(gid, '\t\t%s\n', str);
str = 'double rux = u7;';
fprintf(gid, '\t\t%s\n', str);
str = 'double rvx = u8;';
fprintf(gid, '\t\t%s\n', str);
str = 'double rwx = u9;';
fprintf(gid, '\t\t%s\n', str);
str = 'double rEx = u10;';
fprintf(gid, '\t\t%s\n', str);
str = 'double ry = u11;';
fprintf(gid, '\t\t%s\n', str);
str = 'double ruy = u12;';
fprintf(gid, '\t\t%s\n', str);
str = 'double rvy = u13;';
fprintf(gid, '\t\t%s\n', str);
str = 'double rwy = u14;';
fprintf(gid, '\t\t%s\n', str);
str = 'double rEy = u15;';
fprintf(gid, '\t\t%s\n', str);
str = 'double rz = u16;';
fprintf(gid, '\t\t%s\n', str);
str = 'double ruz = u17;';
fprintf(gid, '\t\t%s\n', str);
str = 'double rvz = u18;';
fprintf(gid, '\t\t%s\n', str);
str = 'double rwz = u19;';
fprintf(gid, '\t\t%s\n', str);
str = 'double rEz = u20;';
fprintf(gid, '\t\t%s\n', str);
fprintf(gid, '\n');  

fid = fopen(temporaryFileName1,'r');    
tline = fgetl(fid); 
i=1; a1 = 0;       
while ischar(tline)        
    str = tline;

    i1 = strfind(str,'[');        
    i2 = strfind(str,']');        
    if isempty(i1)==0    
        a2 = str2num(str((i1+1):(i2-1)));                        
        for j = a1:(a2-1)                
            strj = ['f[' num2str(j) '*ng+i] += 0.0;'];
            fprintf(gid, '\t\t%s\n', strj);                  
        end
        a1 = a2+1;              
    end

    str = strrep(str, '  ', '');
    str = strrep(str, 'A0[', 'f[');
    str = strrep(str, '][0] = ', '*ng+i] = ');   
    str = strrep(str, '] = ', '] += ');   
    if isempty(i1)==1
        str = ['double ' str];
    end

    fprintf(gid, '\t\t%s\n', str);    
    tline = fgetl(fid);        
    i=i+1;   
end
fclose(fid);
fprintf(gid, '\n');  

str = '}';
fprintf(gid, '\t%s\n', str);             

str = 'if (computeJacobian == 1) {';
fprintf(gid, '\n\t%s', str); 

str = 'for (int i = 0; i <ng; i++) {';
fprintf(gid, '\n\t\t%s\n', str);         

% for i=1:length(pg)
%     str = ['double x' num2str(i) ' = pg[' num2str(i-1) '*ng+i];'];
%     fprintf(gid, '\t\t%s\n', str);                  
% end
str = 'double x = pg[0*ng+i];';
fprintf(gid, '\t\t\t%s\n', str);
str = 'double y = pg[1*ng+i];';
fprintf(gid, '\t\t\t%s\n', str);
str = 'double z = pg[2*ng+i];';
fprintf(gid, '\t\t\t%s\n', str);
str = 'double h = pg[3*ng+i];';
fprintf(gid, '\t\t\t%s\n', str);
str = 'double Delta = h;        // TODO: Implement Delta = V^(1/3)';
fprintf(gid, '\t\t\t%s\n', str);
str = 'double C_d = Cs[i];';
fprintf(gid, '\t\t\t%s\n', str);
fprintf(gid, '\t\t\n');
for i=1:length(udg)
    str = ['double u' num2str(i) ' = udg[' num2str(i-1) '*ng+i];'];
    fprintf(gid, '\t\t\t%s\n', str);                  
end
fprintf(gid, '\n');  

str = 'double r = u1;';
fprintf(gid, '\t\t\t%s\n', str);
str = 'double ru = u2;';
fprintf(gid, '\t\t\t%s\n', str);
str = 'double rv = u3;';
fprintf(gid, '\t\t\t%s\n', str);
str = 'double rw = u4;';
fprintf(gid, '\t\t\t%s\n', str);
str = 'double rE = u5;';
fprintf(gid, '\t\t\t%s\n', str);
str = 'double rx = u6;';
fprintf(gid, '\t\t\t%s\n', str);
str = 'double rux = u7;';
fprintf(gid, '\t\t\t%s\n', str);
str = 'double rvx = u8;';
fprintf(gid, '\t\t\t%s\n', str);
str = 'double rwx = u9;';
fprintf(gid, '\t\t\t%s\n', str);
str = 'double rEx = u10;';
fprintf(gid, '\t\t\t%s\n', str);
str = 'double ry = u11;';
fprintf(gid, '\t\t\t%s\n', str);
str = 'double ruy = u12;';
fprintf(gid, '\t\t\t%s\n', str);
str = 'double rvy = u13;';
fprintf(gid, '\t\t\t%s\n', str);
str = 'double rwy = u14;';
fprintf(gid, '\t\t\t%s\n', str);
str = 'double rEy = u15;';
fprintf(gid, '\t\t\t%s\n', str);
str = 'double rz = u16;';
fprintf(gid, '\t\t\t%s\n', str);
str = 'double ruz = u17;';
fprintf(gid, '\t\t\t%s\n', str);
str = 'double rvz = u18;';
fprintf(gid, '\t\t\t%s\n', str);
str = 'double rwz = u19;';
fprintf(gid, '\t\t\t%s\n', str);
str = 'double rEz = u20;';
fprintf(gid, '\t\t\t%s\n', str);
fprintf(gid, '\n');  

fid = fopen(temporaryFileName2,'r');    
tline = fgetl(fid); 
i=1; a1=0;      
while ischar(tline)        
    str = tline;
    %str = strrep(tline, ';', '');
    %str = strrep(str, 't0 = ', '');        

    i1 = strfind(str,'[');        
    i2 = strfind(str,']');        
    if isempty(i1)==0    
        a2 = str2num(str((i1+1):(i2-1)));                        
        for j = a1:(a2-1)                
            strj = ['f_udg[' num2str(j) '*ng+i] += 0.0;'];
            fprintf(gid, '\t\t\t%s\n', strj);  
        end
        a1 = a2+1;              
    end

    str = strrep(str, '  ', '');        
    str = strrep(str, 'A0[', 'f_udg[');
    str = strrep(str, '][0] = ', '*ng+i] = ');  
    str = strrep(str, '] = ', '] += ');
    if isempty(i1)==1
        str = ['double ' str];
    end

    fprintf(gid, '\t\t\t%s\n', str);                          
    tline = fgetl(fid);        
    i=i+1;        
end
fclose(fid);
fprintf(gid, '\n');  

str = '}';
fprintf(gid, '\t\t%s\n', str);             

str = '}';
fprintf(gid, '\t%s\n', str); 

str = '}';
fprintf(gid, '%s\n', str);             
fprintf(gid, '\n');

delete(temporaryFileName1);
delete(temporaryFileName2);
