
fileName = 'flux_ns3dNEW2_';

ncu = 5;
nd = 3;
nc = 20;

syms r ru rv rw rE rx rux rvx rwx rEx ry ruy rvy rwy rEy rz ruz rvz rwz rEz
syms gam epsilon Re Pr Minf tau
syms time
syms mu T p

udg = [r ru rv rw rE rx rux rvx rwx rEx ry ruy rvy rwy rEy rz ruz rvz rwz rEz];
param = [gam epsilon Re Pr Minf tau];

% Field variables:
r1 = 1./r;
u = ru./r;
v = rv./r;
w = rw./r;
E = rE./r;
q = 0.5*(u.*u+v.*v+w.*w);
p = (gam-1)*(rE-r.*q);
H = E+p./r;
T = p./((gam-1)*r);

% Derivatives of field variables:
ux = (rux - u.*rx)./r;
vx = (rvx - v.*rx)./r;
wx = (rwx - w.*rx)./r;
qx  = u.*ux + v.*vx + w.*wx;
px  = (gam-1)*(rEx - rx.*q - r.*qx);
Tx  = (px.*r - p.*rx)./((gam-1)*r.^2);

uy = (ruy - u.*ry)./r;
vy = (rvy - v.*ry)./r;
wy = (rwy - w.*ry)./r;
qy  = u.*uy + v.*vy + w.*wy;
py  = (gam-1)*(rEy - ry.*q - r.*qy);
Ty  = (py.*r - p.*ry)./((gam-1)*r.^2);

uz = (ruz - u.*rz)./r;
vz = (rvz - v.*rz)./r;
wz = (rwz - w.*rz)./r;
qz  = u.*uz + v.*vz + w.*wz;
pz  = (gam-1)*(rEz - rz.*q - r.*qz);
Tz  = (pz.*r - p.*rz)./((gam-1)*r.^2);

% Viscous stresses:
txx = (1/Re)*mu.*(2*ux - (2./3)*(ux+vy+wz));
txy = (1/Re)*mu.*(uy+vx);
txz = (1/Re)*mu.*(uz+wx);
tyx = (1/Re)*mu.*(vx+uy);
tyy = (1/Re)*mu.*(2*vy - (2./3)*(ux+vy+wz));
tyz = (1/Re)*mu.*(vz+wy);
tzx = (1/Re)*mu.*(wx+uz);
tzy = (1/Re)*mu.*(wy+vz);
tzz = (1/Re)*mu.*(2*wz - (2./3)*(ux+vy+wz));

factor1 = (1/Re) * ((gam.*mu)./Pr);

F11 = ru;
F12 = rv;
F13 = rw;
F21 = ru.*u + p + txx;
F22 = ru.*v + txy;
F23 = ru.*w + txz;
F31 = rv.*u + tyx;
F32 = rv.*v + p + tyy;
F33 = rv.*w + tyz;
F41 = rw.*u + tzx;
F42 = rw.*v + tzy;
F43 = rw.*w + p + tzz;
F51 = ru.*H + txx.*u + txy.*v + txz.*w + factor1.*Tx;
F52 = rv.*H + tyx.*u + tyy.*v + tyz.*w + factor1.*Ty;
F53 = rw.*H + tzx.*u + tzy.*v + tzz.*w + factor1.*Tz;

% Fluxes:
F = [F11; F21; F31; F41; F51; F12; F22; F32; F42; F52; F13; F23; F33; F43; F53];

% Jacobian of the fluxes:
F_udg = simplify(jacobian(F,udg));

% F_mu = simplify(jacobian(F,mu));

F = reshape(F,[ncu,nd]);
F_udg = reshape(F_udg,[ncu,nd,nc]);

% Generate C code:
disp('  ')
disp('  ')
disp('Auto-generating C code...'); tic

temporaryFileName1 = ['tmp_', fileName, '1.c'];
temporaryFileName2 = ['tmp_', fileName, '2.c'];

% generate temporary C file
ccode(F(:),'file',temporaryFileName1)
ccode(F_udg(:),'file',temporaryFileName2)

fileNameC = [fileName, '.c'];
delete(fileNameC);
gid = fopen(fileNameC,'wt');

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
            strj = ['f[' num2str(j) '*ng+i] = 0.0;'];
            fprintf(gid, '\t\t%s\n', strj);                  
        end
        a1 = a2+1;              
    end

    str = strrep(str, '  ', '');
    str = strrep(str, 'A0[', 'f[');
    str = strrep(str, '][0] = ', '*ng+i] = ');   
    str = strrep(str, '] = ', '] = ');   
    if isempty(i1)==1
        str = ['double ' str];
    end

    fprintf(gid, '\t\t%s\n', str);    
    tline = fgetl(fid);        
    i=i+1;   
end
fclose(fid);
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
            strj = ['f_udg[' num2str(j) '*ng+i] = 0.0;'];
            fprintf(gid, '\t\t\t%s\n', strj);  
        end
        a1 = a2+1;              
    end

    str = strrep(str, '  ', '');        
    str = strrep(str, 'A0[', 'f_udg[');
    str = strrep(str, '][0] = ', '*ng+i] = ');  
    str = strrep(str, '] = ', '] = ');
    if isempty(i1)==1
        str = ['double ' str];
    end

    fprintf(gid, '\t\t\t%s\n', str);                          
    tline = fgetl(fid);        
    i=i+1;        
end
fclose(fid);

delete(temporaryFileName1);
delete(temporaryFileName2);
