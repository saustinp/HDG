
fileName = 'molecularAV_3d';

ncu = 5;
nd = 3;
nc = 20;

% Input variables:
syms r ru rv rw rE rx rux rvx rwx rEx ry ruy rvy rwy rEy rz ruz rvz rwz rEz
syms ux uy uz vx vy vz wx wy wz av
syms gam epslm Re Minf Pr tau

udg = [r ru rv rw rE rx rux rvx rwx rEx ry ruy rvy rwy rEy rz ruz rvz rwz rEz];
param = [gam epslm Re Pr Minf tau];

% Field variables:
u = ru./r;
v = rv./r;
w = rw./r;
ux = (rux - u.*rx)./r;
vx = (rvx - v.*rx)./r;
wx = (rwx - w.*rx)./r;
uy = (ruy - u.*ry)./r;
vy = (rvy - v.*ry)./r;
wy = (rwy - w.*ry)./r;
uz = (ruz - u.*rz)./r;
vz = (rvz - v.*rz)./r;
wz = (rwz - w.*rz)./r;

txx_AV = r.*av.*(ux + ux - (2.0/3.0)*(ux+vy+wz));
txy_AV = r.*av.*(uy + vx);
txz_AV = r.*av.*(uz + wx);
tyx_AV = r.*av.*(vx + uy);
tyy_AV = r.*av.*(vy + vy - (2.0/3.0)*(ux+vy+wz));
tyz_AV = r.*av.*(vz + wy);
tzx_AV = r.*av.*(wx + uz);
tzy_AV = r.*av.*(wy + vz);
tzz_AV = r.*av.*(wz + wz - (2.0/3.0)*(ux+vy+wz));

% Artificial viscosity fluxes:
F_AV = [0; txx_AV; txy_AV; txz_AV; txx_AV.*u + txy_AV.*v + txz_AV.*w; 0; tyx_AV; tyy_AV; tyz_AV; tyx_AV.*u + tyy_AV.*v + tyz_AV.*w; 0; tzx_AV; tzy_AV; tzz_AV; tzx_AV.*u + tzy_AV.*v + tzz_AV.*w];
% F_AV = [av.*rx; txx_AV; tyx_AV; tzx_AV; txx_AV.*u + txy_AV.*v + txz_AV.*w + av.*rEx; av.*ry; txy_AV; tyy_AV; tzy_AV; tyx_AV.*u + tyy_AV.*v + tyz_AV.*w + av.*rEy; av.*rz; txz_AV; tyz_AV; tzz_AV; tzx_AV.*u + tzy_AV.*v + tzz_AV.*w + av.*rEz];

% Jacobian of the artificial viscosity fluxes:
F_AV_udg = jacobian(F_AV,udg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F_AV = reshape(F_AV,[ncu,nd]);
F_AV_udg = reshape(F_AV_udg,[ncu,nd,nc]);

F = F_AV;
F_udg = F_AV_udg;

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
