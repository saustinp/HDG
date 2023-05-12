function cuda_check(mesh, master, app, SH, UDG, UH, uinf, param, time, fc_q, fc_u, tdep, f2e, nn, nme, nmf)

% cuda_check(mesh, master, app, 0*UDG, UDG, UH, uinf, param, time, 1.0, 1.0, 0, f2e, nn, nme, nmf);

[npe, ncx, ne] = size(mesh.dgnodes);
nge = master.ngv;
nd = master.nd;
nc = size(UDG,2);
ncu = nn(3);
ncq = ncu*nd;
nf = mesh.nf;
ngf = master.ngf;
npf = master.npf;
nfe = size(master.perm,2);

err = zeros(10,1);

un = reshape(permute(UDG,[1 3 2]),[npe ne*nc]);
ug = master.shapvt(:,:,1)*un;

cUDG = zeros(npe*ne, nc);
cUDGg = zeros(nge*ne, nc);
cUDG = cuda_getelemnodes(cUDG, UDG, npe, ne, nc, 1, ne);    
cUDGg = cuda_node2gauss(cUDGg, cUDG, master.shapvt, nge, npe, ne*nc);

err(1) = max(abs(ug(:)-cUDGg(:)));
if max(err(1))<1e-12
    disp('check cuda_node2gauss: success');
else
    disp('check cuda_node2gauss: fail');
end

pn = reshape(permute(mesh.dgnodes,[1 3 2]),[npe ne ncx]);
[pg, Xx, jac] = volgeom(master.shapvt,pn);

cpg = 0*pg;
cXx = 0*Xx;
cjac = 0*jac;
cJg = 0*Xx;
cpn = zeros(npe*ne, ncx);
cpn = cuda_getelemnodes(cpn, mesh.dgnodes, npe, ne, ncx, 1, ne);    
[cpg, cXx, cjac] = cuda_elemgeom(cpg, cXx, cjac, cpn, cJg, master.shapvt, ne, nge, npe, nd, ncx);

err(1) = max(abs(pg(:)-cpg(:)));
err(2) = max(abs(Xx(:)-cXx(:)));
err(3) = max(abs(jac(:)-cjac(:)));
if max(err)<1e-12
    disp('check cuda_elemgeom: success');
else
    disp('check cuda_elemgeom: fail');
end

Cu = zeros(npe*ne,ncu,nd);
fqg = zeros(nge*ne,ncu);
ug = reshape(ug,[nge*ne nc]);
Xx = reshape(Xx,[nge*ne nd nd]);
for i=1:nd
    for m=1:ncu
        fqg(:,m) = ug(:,m).*Xx(:,1,i);
    end
    for j=2:nd
        for m=1:ncu
            fqg(:,m) = fqg(:,m) + ug(:,m).*Xx(:,j,i);
        end
    end    
    Cu(:,:,i) = reshape(master.shapvg(:,:,i+1)*reshape(fqg,[nge ne*ncu]),[npe*ne ncu]);    
end
Cu = permute(reshape(Cu,[npe ne ncq]),[1 3 2]);

shapen = master.shapvt;
shapeg = master.shapvg;
[Rq,tempen,tempeg] = cuda_memalloc(nn);
Rq = cuda_Rq_elem(0*Rq, UDG, mesh.dgnodes, tempen, tempeg, shapen, shapeg, nn, nme);

err(1) = max(abs(Cu(:)-Rq(:)));
if max(err(1))<1e-12
    disp('check cuda_Rq_elem: success');
else
    disp('check cuda_Rq_elem: fail');
end

% shapfn = master.shapft;
% shapfg = master.shapfg;
% n1 = ncx+nd+1+2*nc+ncu+ncu*nd;             
% blkszf = nf;
% tempfg = zeros(ngf*blkszf,n1);
% tempfn = zeros(npf*blkszf,max([ncx,nc]));
% nmf = [1; nf; 0];
% [Rq,cxhg,cnlg,cjcg,ug1,ug2,uhg] = cuda_Rq_face(0*Rq, UDG, mesh.dgnodes, UH, uinf, tempfn, tempfg, shapfn, shapfg, param, time, f2e, nn, nmf);
% cxhg = reshape(cxhg, [ngf nf ncx]);
% cnlg = reshape(cnlg, [ngf nf nd]);
% cjcg = reshape(cjcg, [ngf nf 1]);
% 
% [xhg, nlg, jcg] = facegeom(master.shapft,pn,mesh.perm);
% xhg = reshape(xhg, [ngf nfe ne ncx]);
% nlg = reshape(nlg, [ngf nfe ne nd]);
% jcg = reshape(jcg, [ngf nfe ne 1]);
% 
% for i = 1:nf
%     e = mesh.f(i,end-1);
%     a = mesh.t2f(e,:)==i;        
%     eg = reshape(xhg(:,a,e,:),[ngf ncx])-reshape(cxhg(:,i,:),[ngf ncx]);
%     err(1) = max(err(1),max(abs(eg(:))));
%     eg = reshape(nlg(:,a,e,:),[ngf nd])-reshape(cnlg(:,i,:),[ngf nd]);
%     err(2) = max(err(2),max(abs(eg(:))));
%     eg = reshape(jcg(:,a,e,:),[ngf 1])-reshape(cjcg(:,i,:),[ngf 1]);
%     err(3) = max(err(3),max(abs(eg(:))));
% end
% if max(err)<1e-12
%     disp('check cuda_facegeom: success');
% else
%     disp('check cuda_facegeom: fail');
% end

shapfn = master.shapft;
shapfg = master.shapfg;
n1 = ncx+nd+1+2*nc+ncu+ncu*nd;             
blkszf = nn(12);
tempfg = zeros(ngf*blkszf,n1);
tempfn = zeros(npf*blkszf,max([ncx,nc]));
Rq = cuda_Rq_face(Rq, UDG, mesh.dgnodes, permute(reshape(UH,[ncu npf nf]),[2 1 3]), uinf, tempfn, tempfg, shapfn, shapfg, param, time, f2e, nn, nmf);
Euh = getq2(master, mesh, UDG, UH, SH, fc_q);

err(1) = max(abs(Euh(:)-Rq(:)));
if max(err(1))<1e-12
    disp('check cuda_Rq_face: success');
else
    disp('check cuda_Rq_face: fail');
end

Mi = massinv(master, mesh);
q = cuda_getq(0*Rq, Mi, Rq, SH, UDG, mesh.dgnodes, permute(reshape(UH,[ncu npf nf]),[2 1 3]), uinf, tempen, tempeg, tempfn, tempfg, shapen, shapeg, shapfn, shapfg, param, time, fc_q, f2e, nn, nme, nmf);
QDG = getq(master, mesh, UDG, UH, SH, fc_q);

err(1) = max(abs(q(:)-QDG(:)));
if max(err(1))<1e-10
    disp('check cuda_getq: success');
else
    disp('check cuda_getq: fail');
end

Ru = zeros(npe, ncu, ne);
Ru = cuda_Ru_elem(Ru, UDG, mesh.dgnodes, SH(:,1:ncu,:), tempen, tempeg, shapen, shapeg, param, time, fc_u, tdep, nn, nme);
Rv = volint(master,app,permute(mesh.dgnodes,[1 3 2]),permute(UDG,[1 3 2]),permute(SH,[1 3 2]));

err(1) = max(abs(Ru(:)-Rv(:)));
if max(err(1))<1e-12
    disp('check cuda_Ru_elem: success');
else
    disp('check cuda_Ru_elem: fail');
end

Ru = cuda_Ru_face(Ru, UDG, mesh.dgnodes, uinf(1), tempfn, tempfg, shapfn, shapfg, param, time, f2e, nn, nmf);
app.fbou = 'ldgfbou';
app.fhat = 'ldgfhat';
Rv = faceint(master,mesh,app,UDG,permute(reshape(UH,[ncu npf nf]),[2 1 3]),Rv);

err(1) = max(abs(Ru(:)-Rv(:)));
if max(err(1))<1e-12
    disp('check cuda_Ru_face: success');
else
    disp('check cuda_Ru_face: fail');
end

Ru = cuda_residual(Mi, 0*Ru, 0*q, 0*Rq, UDG, SH, mesh.dgnodes, permute(reshape(UH,[ncu npf nf]),[2 1 3]), uinf(1), tempen, tempeg, tempfn, tempfg, shapen, shapeg, shapfn, shapfg, param, time, fc_q, fc_u, tdep, f2e, nn, nme, nmf);

err(1) = max(abs(Ru(:)-Rv(:)));
if max(err(1))<1e-12
    disp('check cuda_residual: success');
else
    disp('check cuda_residual: fail');
end


return;

function Ru = volint(master,app,dgnodes,UDG,SH)
% VOLINTND compute volume integrals 

% Obtain dgnodes, Jacobian matrix and determinant at Gauss points
[pg, Xx, jac] = volgeom(master.shapmv,dgnodes);

ne   = size(UDG,2);
nd   = master.nd;
npv  = master.npv;
ngv  = master.ngv;

nc   = app.nc;
ncu  = app.ncu;
arg  = app.arg;
time = app.time;
tdep = app.tdep;
fc_u = app.fc_u;
source = str2func(app.source);
flux   = str2func(app.flux);

% Shap functions and derivatives
shapvt = master.shapvt;
shapvg = reshape(master.shapvg,[npv ngv*(nd+1)]);

% DG solution at Gauss points
udgg = reshape(UDG,[npv ne*nc]);
udgg = shapvt(:,:,1)*udgg;
udgg = reshape(udgg,[ngv*ne nc]);

% Fluxes and source at Gauss points
s = source( pg, udgg, arg, time); 
s     = reshape(s(:,1:ncu),[ngv*ne ncu]);
% Update source term for time-dependent problems
if tdep    
    Stn = reshape(SH(:,:,1:ncu),[npv ne*ncu]);
    Stg = shapvt(:,:,1)*Stn;
    Stg = reshape(Stg,[ngv*ne ncu]);

    s = s + Stg - udgg(:,1:ncu)*fc_u;    
end
wrk = zeros(ngv*(nd+1),ne*ncu);
wrk(1:ngv,:) =  reshape(bsxfun(@times,s,jac),[ngv ne*ncu]);

f = flux( pg, udgg, arg, time);
f     = reshape(f,[ngv ne ncu nd]);
for i=1:nd
    fk = bsxfun(@times,f(:,:,:,1),Xx(:,:,1,i));    
    for j=2:nd
        fk = fk + bsxfun(@times,f(:,:,:,j),Xx(:,:,j,i));        
    end
    wrk(i*ngv+1:(i+1)*ngv,:) = reshape(fk,[ngv ne*ncu]);   
end

% Volume residual
% [Phi Phi_xi Phi_eta] x [S.*jac; Fx.*Xx(:,:,1,1)+Fy.*Xx(:,:,2,1); Fx.*Xx(:,:,1,2)+Fy.*Xx(:,:,2,2)]
Ru = shapvg*wrk; % [npv ngv*(nd+1)] x [ngv*(nd+1) ne*ncu] 
Ru = permute(reshape(Ru,[npv ne ncu]),[1 3 2]); 

% tmp = reshape(wrk,[ngv (nd+1) ne ncu]);
% squeeze(tmp(:,:,1,1))

function RU = faceint(master,mesh,app,UDG,UH,RU)
% FACEINTND compute face integrals 

npf = size(mesh.perm,1);
fbou   = str2func(app.fbou);
fhat   = str2func(app.fhat);

% Shap functions 
perm            = master.perm(:,:,1);
shapft          = master.shapft(:,:,1);
shapfg          = master.shapfg(:,:,1);
elcon = reshape(mesh.elcon,[npf mesh.nfe mesh.ne]);
for i = 1:mesh.nf
    fi = mesh.f(i,end-1:end); % obtain two elements sharing the same face i      
    if fi(2)>0           % face i is an interior face                
        kf = mesh.t2f(fi,:);         % obtain neighboring faces 
        i1 = kf(1,:)==i;  % obtain the index of face i in the 1st element
        i2 = kf(2,:)==i;  % obtain the index of face i in the 2nd element                                            
        j1 = elcon(:,i1,fi(1)) - (i-1)*npf;        
        j2 = elcon(:,i2,fi(2)) - (i-1)*npf;                
        udg1 = shapft*UDG(perm(j1,i1),:,fi(1));                
        udg2 = shapft*UDG(perm(j2,i2),:,fi(2));                
        uhg  = shapft*UH(:,:,i);
        xdg1 = mesh.dgnodes(perm(j1,i1),:,fi(1));
        [pg1, nlg1, jac1] = facegeom(master.shapmf,xdg1);
        
        %FH = fhat(nlg, pg, udgg, uhg, arg, time);   
        fhg = fhat(nlg1, pg1, udg1, udg2, uhg, app.arg, app.time);
        cnt = shapfg*diag(jac1)*fhg;                
        RU(perm(j1,i1),:,fi(1)) = RU(perm(j1,i1),:,fi(1)) - cnt;
        RU(perm(j2,i2),:,fi(2)) = RU(perm(j2,i2),:,fi(2)) + cnt;          
    else % face i is a boundary face
        kf = mesh.t2f(fi(1),:); % obtain neighboring faces 
        i1 = kf(1,:)==i;  % obtain the index of face i in the 1st element                 
        j1 = elcon(:,i1,fi(1)) - (i-1)*npf;        
        udg1 = shapft*UDG(perm(j1,i1),:,fi(1));
        uhg  = shapft*UH(:,:,i);
        xdg1 = mesh.dgnodes(perm(j1,i1),:,fi(1));
        [pg1, nlg1, jac1] = facegeom(master.shapmf,xdg1);                        
        b = -fi(2);
        ib = app.bcm(b);
        uinf = app.bcs(b,:);
        fhg = fbou(ib, uinf, nlg1, pg1, udg1, uhg, app.arg, app.time);
        cnt = shapfg*diag(jac1)*fhg;               
        RU(perm(j1,i1),:,fi(1)) = RU(perm(j1,i1),:,fi(1)) - cnt;        
    end        
end 

function [pg, nlg, jac] = facegeom(shapgeomft,pn)
% FACEGEOM computes dg nodes, Jacobian determinant and normal vectors at Gauss points  

%   [pg, nlg, jac] = facegeom(shapgeomft,dgnodes,perm)
%
%    SHAPGEOMFT :  Shape functions and derivatives at Gauss points
%    DGNODES    :  Geometry DG nodes 
%    PERM       :  Indices of the boundary nodes
%    PG         :  Physical nodes at Gauss points 
%    nlg        :  Normal vector at Gauss points
%    jac        :  Determinant of the Jacobian mapping 

nq    = size(pn,2);
ngf   = size(shapgeomft,1);
npf   = size(shapgeomft,2);
nd    = size(shapgeomft,3);

if nd>1
    dshapft  = reshape(permute(shapgeomft(:,:,2:nd),[1 3 2]),[ngf*(nd-1) npf]);
    dpg = dshapft*pn(:,1:nd);
    dpg = permute(reshape(dpg,[ngf nd-1 nd]), [1 3 4 5 2]);    
    dpg = reshape(dpg,[ngf,nd,nd-1]);    
end

shapgeomft   = shapgeomft(:,:,1);
pg = shapgeomft*pn;
pg = reshape(pg,[ngf nq]);

switch nd
    case 1
        jac = ones(1,1);
        nlg = [ones(1,1); -ones(1,1)];        
        nlg = nlg(:);
    case 2
        jac = sqrt(dpg(:,1).^2+dpg(:,2).^2);
        nlg   = [dpg(:,2),-dpg(:,1)];
        nlg   = bsxfun(@rdivide, nlg, jac);
    case 3
        nlg(:,1) = dpg(:,2,1).*dpg(:,3,2) - dpg(:,3,1).*dpg(:,2,2);
        nlg(:,2) = dpg(:,3,1).*dpg(:,1,2) - dpg(:,1,1).*dpg(:,3,2);
        nlg(:,3) = dpg(:,1,1).*dpg(:,2,2) - dpg(:,2,1).*dpg(:,1,2);
        jac = sqrt(nlg(:,1).^2+nlg(:,2).^2+nlg(:,3).^2);
        nlg   = bsxfun(@rdivide, nlg, jac);
    otherwise
        error('Dimension is not implemented');
end

function QDG = getq2(master, mesh, UDG, UH, SH, fc_q)

% MiCE: Cell with Mi, C and E for all elements. If used as output variable,
% very large memory will be required.

if nargin < 5; SH = []; end
if nargin < 6; fc_q = 1; end

% get dimensions
nd    = master.nd;
npv   = master.npv;
ngv   = master.ngv;
ngf   = master.ngf;
npf   = master.npf;
nfe   = size(master.perm,2);
nch   = size(UH,1);
ncq   = nch*nd;
ne    = size(UDG,3);

%UH    = UH(:,mesh.elcon);
UH = permute(reshape(UH(:,mesh.elcon),[nch nfe*npf ne]),[2 1 3]);

% get shape functions and their derivatives
perm = master.perm;
dshapvt = reshape(permute(master.shapvt(:,:,2:nd+1),[1 3 2]),[ngv*nd npv]);
shapvgdotshapvl = reshape(master.shapvgdotshapvl,[npv*npv ngv (nd+1)]);
dshapft = reshape(permute(master.shapft(:,:,2:nd),[1 3 2]),[ngf*(nd-1) npf]);
shapfgdotshapfc = master.shapfgdotshapfc(:,:,1);
shapnvl = mkshape(mesh.porder,master.plocvl,master.plocvl,mesh.elemtype);
dshapnvt = reshape(permute(shapnvl(:,:,2:nd+1),[2 3 1]),[npv*nd npv]);

% allocate memory
QDG = zeros(npv,ncq,ne);

% loop through each element        
for k=1:1:ne 
%    if rem(k-1,10000) == 0; disp(['Element No. ', num2str(k),' / ', num2str(ne)]); end

    % get dg nodes
    dg = mesh.dgnodes(:,1:nd,k);    

    % get dg nodes on faces
    pn = reshape(dg(perm,:,:),[npf nfe*nd]);

    % compute volumetic Jacobian matrix 
    Jg = dshapvt*dg(:,1:nd);
    Jg = reshape(Jg,[ngv nd nd]);  
    dshapnvt = reshape(dshapnvt,[npv*nd, npv]);
    % compute the determintant and inverse
    [jac,Xx] = volgeom2(Jg);

    dpg = dshapft*pn;
    dpg = permute(reshape(dpg,[ngf nd-1 nfe nd]), [1 3 4 2]);    
    dpg = reshape(dpg,[ngf*nfe,nd,nd-1]); 
    % get the normal vectors and Jacobian determinants on faces
    [nlg, jacf] = facegeom2(dpg,nd);

    % mass matrix
    M = reshape(shapvgdotshapvl(:,:,1)*jac,[npv npv]);    
    Mi= inv(M)/fc_q;

    % convection matrices
    C = zeros(npv,npv,nd); 
    for i=1:nd
        C(:,:,i) = reshape(shapvgdotshapvl(:,:,2)*Xx(:,i,1),[npv npv]);        
        for j=2:nd
            C(:,:,i) = C(:,:,i) + reshape(shapvgdotshapvl(:,:,j+1)*Xx(:,i,j),[npv npv]);
        end            
    end    
    C = reshape(permute(C,[1 3 2]),[npv*nd npv]);

    % face matrices
    E = zeros(npv,npf*nfe,nd);    
    for i=1:nd    
        njc = reshape(nlg(:,i).*(-jacf),[ngf,nfe]);
        wrk = reshape(shapfgdotshapfc*njc,[npf npf*nfe]);   
        for j=1:nfe
            E(perm(:,j),(j-1)*npf+1:j*npf,i) = wrk(1:npf,(j-1)*npf+1:j*npf,:);      
        end            
    end    
    E = reshape(permute(E,[1 3 2]),[npv*nd npf*nfe]);

    % obtain the current solution and numerical trace    
    u    = UDG(:,1:nch,k);    
    uh   = UH(:,:,k);     
    %sh   = SH(:,end-ncq+1:end,k);        
    
    % compute dq
    Cu   = reshape(C*u,[npv nd*nch]);
    Euh  = reshape(E*uh,[npv nd*nch]);
    MiCu = permute(reshape(Mi*Cu,[npv nd nch]), [1 3 2]);
    MiEu = permute(reshape(Mi*Euh,[npv nd nch]), [1 3 2]);
    q    = reshape(MiEu,[npv nch*nd]) - reshape(MiCu,[npv nch*nd]);         
    q = Euh - Cu;
    QDG(:,:,k) = q; 
end

if ~isempty(SH) && fc_q ~= 0
    QDG = QDG + SH(:,end-ncq+1:end,:)/fc_q;
end

function [jac,Xx] = volgeom2(Jg)

ngv = size(Jg,1); 
nd  = size(Jg,2);
switch nd
    case 1
        jac = Jg;
        Xx = -ones(ngv,1);
    case 2
        jac = Jg(:,1,1).*Jg(:,2,2) - Jg(:,1,2).*Jg(:,2,1);
        Xx(:,1,1) = -Jg(:,2,2);
        Xx(:,2,1) = Jg(:,2,1);
        Xx(:,1,2) = Jg(:,1,2);
        Xx(:,2,2) = -Jg(:,1,1);
    case 3
        jac = Jg(:,1,1).*Jg(:,2,2).*Jg(:,3,3) - Jg(:,1,1).*Jg(:,3,2).*Jg(:,2,3)+ ...
              Jg(:,2,1).*Jg(:,3,2).*Jg(:,1,3) - Jg(:,2,1).*Jg(:,1,2).*Jg(:,3,3)+ ...
              Jg(:,3,1).*Jg(:,1,2).*Jg(:,2,3) - Jg(:,3,1).*Jg(:,2,2).*Jg(:,1,3);            
        Xx(:,1,1) = Jg(:,2,3).*Jg(:,3,2) - Jg(:,2,2).*Jg(:,3,3);
        Xx(:,2,1) = Jg(:,2,1).*Jg(:,3,3) - Jg(:,2,3).*Jg(:,3,1);
        Xx(:,3,1) = Jg(:,2,2).*Jg(:,3,1) - Jg(:,2,1).*Jg(:,3,2);
        Xx(:,1,2) = Jg(:,1,2).*Jg(:,3,3) - Jg(:,1,3).*Jg(:,3,2);
        Xx(:,2,2) = Jg(:,1,3).*Jg(:,3,1) - Jg(:,1,1).*Jg(:,3,3);
        Xx(:,3,2) = Jg(:,1,1).*Jg(:,3,2) - Jg(:,1,2).*Jg(:,3,1);
        Xx(:,1,3) = Jg(:,1,3).*Jg(:,2,2) - Jg(:,1,2).*Jg(:,2,3);
        Xx(:,2,3) = Jg(:,1,1).*Jg(:,2,3) - Jg(:,1,3).*Jg(:,2,1);
        Xx(:,3,3) = Jg(:,1,2).*Jg(:,2,1) - Jg(:,1,1).*Jg(:,2,2);
    otherwise
        error('Dimension is not implemented');
end

function [nlg, jac] = facegeom2(dpg,nd)

switch nd
    case 1
        jac = 1;
        nlg = [1; -1];        
        nlg = nlg(:);
    case 2
        jac = sqrt(dpg(:,1).^2+dpg(:,2).^2);
        nlg   = [dpg(:,2),-dpg(:,1)];
        nlg   = bsxfun(@rdivide, nlg, jac);
    case 3
        nlg(:,1) = dpg(:,2,1).*dpg(:,3,2) - dpg(:,3,1).*dpg(:,2,2);
        nlg(:,2) = dpg(:,3,1).*dpg(:,1,2) - dpg(:,1,1).*dpg(:,3,2);
        nlg(:,3) = dpg(:,1,1).*dpg(:,2,2) - dpg(:,2,1).*dpg(:,1,2);
        jac = sqrt(nlg(:,1).^2+nlg(:,2).^2+nlg(:,3).^2);
        nlg   = bsxfun(@rdivide, nlg, jac);
    otherwise
        error('Dimension is not implemented');
end


