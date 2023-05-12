syms u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12 u13 u14 u15 u16 u17 u18 u19 u20
syms uh1 uh2 uh3 uh4 uh5
syms uinf1 uinf2 uinf3 uinf4 uinf5
syms x1 x2 x3 
syms nl1 nl2 nl3
syms time
syms param1 param2 param3 param4 param5 param6

includegetan = 1;
appname = 'ns';

param = [param1 param2 param3 param4 param5 param6];
if nd==2    
    udg = [u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12];
    uh = [uh1 uh2 uh3 uh4];
    uinf = [uinf1 uinf2 uinf3 uinf4];
    pg = [x1 x2];  
    nl = [nl1 nl2];    
else
    udg = [u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12 u13 u14 u15 u16 u17 u18 u19 u20];
    uh = [uh1 uh2 uh3 uh4 uh5];
    uinf = [uinf1 uinf2 uinf3 uinf4 uinf5];
    pg = [x1 x2 x3];  
    nl = [nl1 nl2 nl3];    
end

gam  = param(1);
gam1 = param(1) - 1.0;          
Re   = param(3);
Pr   = param(4);
Minf = param(5);
tau  = param(6);
Re1  = 1/Re;
M2   = Minf^2;
c23  = 2.0/3.0;
fc   = 1/(gam1*M2*Re*Pr);

ncu = length(uh);
nc = length(udg);
nch = ncu;

if (nd==2) && (axissym==1) % NS in axisymetric coordinates                                              
    r    = udg(1);
    ru   = udg(2);
    rv   = udg(3);
    rE   = udg(4);

    rx   = udg(5);
    rux  = udg(6);
    rvx  = udg(7);
    rEx  = udg(8);

    ry   = udg(9);
    ruy  = udg(10);
    rvy  = udg(11);
    rEy  = udg(12);
    
    xz = pg(1);
    xr = pg(2);

    r1   = 1/r;
    uv   = ru*r1;
    vv   = rv*r1;
    E    = rE*r1;
    q   = 0.5*(uv*uv+vv*vv);
    p    = gam1*(rE-r*q);
    h    = E+p*r1;
    
    fi   = [ru, ru*uv+p, rv*uv,   ru*h, ...
            rv, ru*vv,   rv*vv+p, rv*h];
     
    uv_r  = -uv*r1;
    uv_ru =  r1;
    vv_r  = -vv*r1;
    vv_rv =  r1;

    ux  = (rux - rx*uv)*r1;
    vx  = (rvx - rx*vv)*r1;    
    qx  = uv*ux + vv*vx;
    px  = gam1*(rEx - rx*q - r*qx);
    Tx  = gam*M2*(px*r - p*rx)*r1^2;

    uy  = (ruy - ry*uv)*r1;
    vy  = (rvy - ry*vv)*r1;    
    qy  = uv*uy + vv*vy;
    py  = gam1*(rEy - ry*q - r*qy);
    Ty  = gam*M2*(py*r - p*ry)*r1^2;

%     txx = Re1*c23*(2*ux - vy);
%     txy = Re1*(uy + vx);
%     tyy = Re1*c23*(2*vy - ux);
    
    txx = Re1*c23*(2*ux - vy + vv/xr); % z z
    txy = Re1*(uy + vx);               % z r    
    tyy = Re1*c23*(2*vy - ux + vv/xr); % r r        
    tzz = Re1*c23*(-2*vv/xr - ux - vy); % theta theta
        
    fv = [0, txx, txy, uv*txx + vv*txy + fc*Tx, ...
          0, txy, tyy, uv*txy + vv*tyy + fc*Ty];

    f = fi+fv;
    
    s = -[rv, ru*vv+txy,  rv*vv+tyy-tzz, rv*h + uv*txy + vv*tyy + fc*Ty];
    s = s./xr;      
elseif (nd==2) && (axissym==0) % NS in polar coordinates                                                  
    r    = udg(1);
    ru   = udg(2);
    rv   = udg(3);
    rE   = udg(4);

    rx   = udg(5);
    rux  = udg(6);
    rvx  = udg(7);
    rEx  = udg(8);

    ry   = udg(9);
    ruy  = udg(10);
    rvy  = udg(11);
    rEy  = udg(12);
    
    xr = pg(1); % r
    xt = pg(2); % theta
    
    r1   = 1/r;
    uv   = ru*r1;
    vv   = rv*r1;
    E    = rE*r1;
    q   = 0.5*(uv*uv+vv*vv);
    p    = gam1*(rE-r*q);
    h    = E+p*r1;
    
    fi   = [ru, ru*uv+p, rv*uv,   ru*h, ...
            rv, ru*vv,   rv*vv+p, rv*h];
     
    uv_r  = -uv*r1;
    uv_ru =  r1;
    vv_r  = -vv*r1;
    vv_rv =  r1;

    ux  = (rux - rx*uv)*r1;
    vx  = (rvx - rx*vv)*r1;    
    qx  = uv*ux + vv*vx;
    px  = gam1*(rEx - rx*q - r*qx);
    Tx  = gam*M2*(px*r - p*rx)*r1^2;

    uy  = (ruy - ry*uv)*r1;
    vy  = (rvy - ry*vv)*r1;    
    qy  = uv*uy + vv*vy;
    py  = gam1*(rEy - ry*q - r*qy);
    Ty  = gam*M2*(py*r - p*ry)*r1^2;

    txx = Re1*c23*(2*ux - (vy - uv)/xr); % r r
    txy = Re1*((uy+vv)/xr + vx);  % r theta    
    tyy = Re1*c23*(2*(vy-uv)/xr - ux); % theta theta
        
    fv = [0, txx, txy, uv*txx + vv*txy + fc*Tx, ...
          0, txy, tyy, uv*txy + vv*tyy + fc*Ty/xr];

    f = fi+fv;
    f(5:8) = f(5:8)./xr;
    
    s = -[ru, ru*uv-rv*vv+txx-tyy, 3*rv*uv+3*txy, ru*h + uv*txx + vv*txy + fc*Tx];
    s = s./xr;          
else  % NS in cylindrical coordinates                                                                                     
    r    = udg(1);
    ru   = udg(2); % z
    rv   = udg(3); % r
    rw   = udg(4); % theta
    rE   = udg(5);

    rx   = udg(6);
    rux  = udg(7);
    rvx  = udg(8);
    rwx  = udg(9);
    rEx  = udg(10);

    ry   = udg(11);
    ruy  = udg(12);
    rvy  = udg(13);
    rwy  = udg(14);
    rEy  = udg(15);
    
    rz   = udg(16);
    ruz  = udg(17);
    rvz  = udg(18);
    rwz  = udg(19);
    rEz  = udg(20);
    
    xz = pg(1);
    xr = pg(2);
    xt = pg(3);
    
    r1   = 1/r;
    uv   = ru*r1;
    vv   = rv*r1;
    wv   = rw*r1;
    E    = rE*r1;
    q   = 0.5*(uv*uv+vv*vv+wv*wv);
    p    = gam1*(rE-r*q);
    h    = E+p*r1;
    
    fi = [ru, ru*uv+p, rv*uv,   rw*uv,   ru*h,  ...
          rv, ru*vv,   rv*vv+p, rw*vv,   rv*h, ...
          rw, ru*wv,   rv*wv,   rw*wv+p, rw*h];    
      
    uv_r  = -uv*r1;
    uv_ru =  r1;
    vv_r  = -vv*r1;
    vv_rv =  r1;
    wv_r  = -wv*r1;
    wv_rv =  r1;

    % p = rho*T/(gamma*Minf*Minf)
    % T = gamma*Minf*Minf*p/rho
    % z
    ux  = (rux - rx*uv)*r1;
    vx  = (rvx - rx*vv)*r1;
    wx  = (rwx - rx*wv)*r1;
    qx  = uv*ux + vv*vx + wv*wx;
    px  = gam1*(rEx - rx*q - r*qx);
    Tx  = gam*M2*(px*r - p*rx)*r1^2;

    % r
    uy  = (ruy - ry*uv)*r1;
    vy  = (rvy - ry*vv)*r1;
    wy  = (rwy - ry*wv)*r1;
    qy  = uv*uy + vv*vy + wv*wy;
    py  = gam1*(rEy - ry*q - r*qy);
    Ty  = gam*M2*(py*r - p*ry)*r1^2;

    % theta
    uz  = (ruz - rz*uv)*r1;
    vz  = (rvz - rz*vv)*r1;
    wz  = (rwz - rz*wv)*r1;
    qz  = uv*uz + vv*vz + wv*wz;
    pz  = gam1*(rEz - rz*q - r*qz);
    Tz  = gam*M2*(pz*r - p*rz)*r1^2;
    
    txx = Re1*c23*(2*ux - vy - (wz-vv)/xr); % z z
    txy = Re1*(uy + vx); % z r
    txz = Re1*(uz/xr + wx); % z theta
    tyy = Re1*c23*(2*vy - ux - (wz-vv)/xr); % r r    
    tyz = Re1*((vz+wv)/xr + wy); % r theta        
    tzz = Re1*c23*(2*(wz-vv)/xr - ux - vy); % theta theta
    
    fv = [0, txx, txy, txz, uv*txx + vv*txy + wv*txz + fc*Tx, ...
          0, txy, tyy, tyz, uv*txy + vv*tyy + wv*tyz + fc*Ty,...
          0, txz, tyz, tzz, uv*txz + vv*tyz + wv*tzz + fc*Tz/xr];
       
    f = fi+fv;         
    f(11:15) = f(11:15)./xr;
    
    s = -[rv, ru*vv+txy,  rv*vv-rw*wv+tyy-tzz, 2*rw*vv+2*tyz, rv*h + uv*txy + vv*tyy + wv*tyz + fc*Ty];
    s = s./xr;    
end

if nd==3
    filename1 = ['flux_' appname num2str(nd) 'd' '.m'];
    filename2 = ['source_' appname num2str(nd) 'd'  '.m'];
elseif axissym==1
    filename1 = ['flux_axis' appname num2str(nd) 'd' '.m'];
    filename2 = ['source_axis' appname num2str(nd) 'd'  '.m'];
else
    filename1 = ['flux_polar' appname num2str(nd) 'd' '.m'];
    filename2 = ['source_polar' appname num2str(nd) 'd'  '.m'];
end
    
genmatlabcode;

%filename3 = ['fhat_' appname num2str(nd) 'd' '.c'];

% genccode; % generate source codes
% 
% % generate an application file
% gid = fopen('fluxes.c','w');
% if includegetan==1
%     str = '#include "../getan.c"';
%     fprintf(gid, '%s\n', str);    
% end
% fid = fopen('../../fluxes_template.c','r');    
% tline = fgetl(fid); 
% while ischar(tline)        
%     str = tline;
%     str = strrep(str, 'appname', appname);    
%     fprintf(gid, '%s\n', str);    
%     tline = fgetl(fid);            
% end
% fclose(fid);
% fclose(gid);
% 
% 


