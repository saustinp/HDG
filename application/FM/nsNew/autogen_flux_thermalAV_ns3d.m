
fileName = 'thermalAV_3d';

ncu = 5;
nd = 3;
nc = 20;

% Input variables:
syms r ru rv rw rE rx rux rvx rwx rEx ry ruy rvy rwy rEy rz ruz rvz rwz rEz
syms ux uy uz vx vy vz wx wy wz av
syms gam gam1 epslm Re Minf Pr tau

udg = [r ru rv rw rE rx rux rvx rwx rEx ry ruy rvy rwy rEy rz ruz rvz rwz rEz];
param = [gam epslm Re Pr Minf tau];

% Field variables:
u = ru./r;
v = rv./r;
w = rw./r;
q = 0.5*(u.*u+v.*v+w.*w);
p = gam1*(rE-r.*q);

% Derivatives of field variables:
ux = (rux - u.*rx)./r;
vx = (rvx - v.*rx)./r;
wx = (rwx - w.*rx)./r;
qx  = u.*ux + v.*vx + w.*wx;
px  = gam1*(rEx - rx.*q - r.*qx);
Tx  = (px.*r - p.*rx)./(gam1*r.^2);

uy = (ruy - u.*ry)./r;
vy = (rvy - v.*ry)./r;
wy = (rwy - w.*ry)./r;
qy  = u.*uy + v.*vy + w.*wy;
py  = gam1*(rEy - ry.*q - r.*qy);
Ty  = (py.*r - p.*ry)./(gam1*r.^2);

uz = (ruz - u.*rz)./r;
vz = (rvz - v.*rz)./r;
wz = (rwz - w.*rz)./r;
qz  = u.*uz + v.*vz + w.*wz;
pz  = gam1*(rEz - rz.*q - r.*qz);
Tz  = (pz.*r - p.*rz)./(gam1*r.^2);

% Artificial viscosity fluxes:
F_AV = [0; 0; 0; 0; r.*av.*Tx; 0; 0; 0; 0; r.*av.*Ty; 0; 0; 0; 0; r.*av.*Tz];

% Jacobian of the artificial viscosity fluxes:
F_AV_udg = jacobian(F_AV,udg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate C code:

F_AV = reshape(F_AV,[ncu,nd]);
F_AV_udg = reshape(F_AV_udg,[ncu,nd,nc]);

F = F_AV;
F_udg = F_AV_udg;

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

delete(temporaryFileName1);
delete(temporaryFileName2);
