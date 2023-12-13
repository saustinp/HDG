
function [Ru, Rq, Rh, BD, F, E, G, H, Ju, Jq, Jh, wrb] = faceint(master,app,bf,dgnodes,UDG,UH,Ru,Rq,BD,Ju,Jq)
% FACEINTND compute face integrals 

% Obtain dgnodes, Jacobian matrix and determinant at Gauss points
[pg, nlg, jac] = facegeom(master.shapmf,dgnodes,master.permgeom);

ne   = size(UDG,2);
nd   = master.nd;
npv  = master.npv;
npf  = master.npf;
ngf  = master.ngf;
nfe  = size(master.perm,2);
nc   = app.nc;
ncu  = app.ncu;
nch  = app.nch;
arg  = app.arg;
time = app.time;
bcm  = app.bcm;        
bcs  = app.bcs;
bcd  = app.bcd;        
bcv  = app.bcv;
localsolve = app.localsolve;
adjoint = app.adjoint;
fbou   = str2func(app.fbou);
fhat   = str2func(app.fhat);

% Shap functions 
perm            = master.perm(:,:,1);
shapft          = master.shapft(:,:,1);
shapfg          = master.shapfg(:,:,1);
shapfgdotshapfc = master.shapfgdotshapfc(:,:,1);

if localsolve==0
    E = [];
else    
    tmp = zeros(npv,npf*nfe,ne);
    E = zeros(npv,npf*nfe,nd,ne);        
    for i=1:nd    
        njc = reshape(nlg(:,i).*(jac),[ngf,nfe*ne]);
        wrk = reshape(shapfgdotshapfc*njc,[npf npf*nfe ne]);   
        for j=1:nfe
            tmp(perm(:,j),(j-1)*npf+1:j*npf,:) = wrk(1:npf,(j-1)*npf+1:j*npf,:);      
        end    
        for k=1:ne
            E(:,:,i,k) = tmp(:,:,k);        
        end
    end
    for i = 1:nd
        Rq(:,:,i,:) = Rq(:,:,i,:) - reshape(mapContractK(E(:,:,i,:),UH,1,[2 3],4,1,3,2),[npv ncu 1 ne]);
    end
end
    
% DG solution at Gauss points
udgn = reshape(UDG(perm,:,:),[npf nfe*ne*nc]);
udgg = shapft*udgn;
udgg = reshape(udgg,[ngf*nfe*ne nc]);

uh  = reshape(UH,[npf nfe*ne*nch]);
uhg = shapft*uh;
uhg = reshape(uhg,[ngf*nfe*ne nch]);

% tmp = reshape(uhg,[ngf*nfe ne nch]);
% squeeze(tmp(:,1,:))
% tmp = reshape(udgg,[ngf*nfe ne nc]);
% squeeze(tmp(:,1,:))
% pause

% arg{end} = 0;
% uhg = 0*uhg;
% udgg = 0*udgg;

% Compute numerical fluxes 
[FH, FH_udg, FH_uh] = fhat(nlg, pg, udgg, uhg, arg, time);

% Code to compare C code (DIGASO) and Matlab implementation of fhat
% Flag implemented so we don't have to uncomment it every time
if isfield(app,'check_fhat_flg')
    nch=ncu;
    ncd=3;
    [ng,nc] = size(udgg);
    arg_dig = cell2mat(arg);

    fh_dig = zeros(size(FH));
    fh_udg_dig = zeros(size(FH_udg));
    fh_uh_dig = zeros(size(FH_uh));

    % Flatten all arrays
    nl_dig = reshape(nlg,1,[]);
    pg_dig = reshape(pg,1,[]);
    udgg_dig = reshape(udgg,1,[]);
    uhg_dig = reshape(uhg,1,[]);
    fh_dig = reshape(fh_dig,1,[]);
    fh_udg_dig = reshape(fh_udg_dig,1,[]);
    fh_uh_dig = reshape(fh_uh_dig,1,[]);

    a1=zeros(1,size(fh_dig,2));
    a2=zeros(1,size(fh_udg_dig,2));
    a3=zeros(1,size(fh_uh_dig,2));
    a4=zeros(1,size(pg_dig,2));
    a5=zeros(1,size(udgg_dig,2));
    a6=zeros(1,size(uhg_dig,2));
    a7=zeros(1,size(nl_dig,2));
    a8=zeros(1,size(arg_dig,2));

    disp('Generating MEX file for fhat ...')
    codegen fhat_digaso_matlab -args {a1,a2,a3,a4,a5,a6,a7,a8,1,1,1,1,1,1} ../digaso_compare/fhat.c

    [fh_dig,fh_udg_dig, fh_uh_dig] = fhat_digaso_matlab_mex(fh_dig, fh_udg_dig, fh_uh_dig, pg_dig, udgg_dig, uhg_dig, nl_dig, arg_dig, time, ng, nc, nch, nd, ncd);

    % Reshape all arrays
    fh_dig = reshape(fh_dig, size(FH));
    fh_udg_dig = reshape(fh_udg_dig, size(FH_udg));
    fh_uh_dig = reshape(fh_uh_dig, size(FH_uh));

    str=sprintf('Linf norm in fhat: %d', max(max(fh_dig-FH)));
    disp(str)
    str=sprintf('Linf norm in fhat_udg: %d', max(max(max(fh_udg_dig-FH_udg))));
    disp(str)
    str=sprintf('Linf norm in fhat_uh: %d', max(max(max(fh_uh_dig - FH_uh))));
    disp(str)

    error('Done checking fhat...')
end

% squeeze(pg)
% squeeze(nlg)
% squeeze(udgg)
% squeeze(uhg)
% squeeze(FH)
% pause

% compute boundary fluxes
% fbou   = str2func(app.fbou);
% tm     = find(bf<0);
% nbf    = length(tm(:));
% JH    = zeros(ngf*nbf,ncu);
% JUDG   = zeros(ngf*nbf,nc);
% BH     = zeros(ngf*nbf,nch);
% BH_udg = zeros(ngf*nbf,nch,nc);
% BH_uh  = zeros(ngf*nbf,nch,nch);
% fbou   = str2func(app.fbou);
% an     = zeros(nbf*npf,1);
% bn     = zeros(nbf*ngf,1);
% j      = 1; 
% for k=1:ne                                
%     for is=1:nfe                            
%         if bf(is,k) < 0 % for each boundary face                        
%             im = (j-1)*ngf+1:j*ngf;            
%             in = (k-1)*nfe*ngf+(is-1)*ngf+1:(k-1)*nfe*ngf+is*ngf;            
%             bn(im) = in; 
%             ib = bcm(-bf(is,k));        
%             bv = repmat(bcs(-bf(is,k),:),[ngf 1]);            
%             [fh, dfh_dudg, dfh_duh] = fbou(ib,bv,nlg(in,:),pg(in,:),udgg(in,:),uhg(in,:),arg,time);                                                 
%             ia = bcd(-bf(is,k));        
%             bu = repmat(bcv(-bf(is,k),:),[ngf 1]);            
%             [~,djh_dudg,djh_duh] = fbou_adjoint(ia,bu,nlg(in,:),pg(in,:),udgg(in,:),uhg(in,:),arg,time);                           
%             JUDG(im,:) = djh_dudg;
%             JH(im,:) = djh_duh;
%             BH(im,:) = fh;
%             BH_udg(im,:,:) = dfh_dudg;
%             BH_uh(im,:,:) = dfh_duh;              
%             im = (j-1)*npf+1:j*npf;            
%             in = (k-1)*nfe*npf+(is-1)*npf+1:(k-1)*nfe*npf+is*npf;            
%             an(im) = in; 
%             j = j+1;
%         end                
%     end
% end

tm     = find(bf<0);
nbf    = length(tm(:));
BH     = zeros(ngf*nbf,nch);
BH_udg = zeros(ngf*nbf,nch,nc);
BH_uh  = zeros(ngf*nbf,nch,nch);
JH     = zeros(ngf*nbf,ncu);
JUDG   = zeros(ngf*nbf,nc);
an     = zeros(nbf*npf,1);
bn     = zeros(nbf*ngf,1);
j = 1;
for i = 1:length(bcm) 
    [I,J] = find(bf==-i);
    nfb = length(I);
    pgb  = zeros(ngf*nfb, size(pg,2));
    nlgb = zeros(ngf*nfb, nd);
    udgb = zeros(ngf*nfb, nc);
    uhgb = zeros(ngf*nfb, ncu);
    im = zeros(ngf*nfb,1);
    for k = 1:nfb
        i0 = (j-1)*ngf+1:j*ngf;
        i1 = (k-1)*ngf+1:k*ngf; 
        i2 = (J(k)-1)*nfe*ngf+(I(k)-1)*ngf+1:(J(k)-1)*nfe*ngf+I(k)*ngf; 
        bn(i0) = i2;
        im(i1)  = i0;
        pgb(i1,:)  = pg(i2,:);
        nlgb(i1,:) = nlg(i2,:);
        udgb(i1,:) = udgg(i2,:);
        uhgb(i1,:) = uhg(i2,:);
        i0 = (j-1)*npf+1:j*npf;
        i2 = (J(k)-1)*nfe*npf+(I(k)-1)*npf+1:(J(k)-1)*nfe*npf+I(k)*npf; 
        an(i0) = i2;
        j = j + 1;
    end  
    
    ib = bcm(i);
    bv = repmat(bcs(i,:),[ngf*nfb 1]);
    [fh, dfh_dudg, dfh_duh] = fbou(ib,bv,nlgb,pgb,udgb,uhgb,arg,time);
    
    % Code to compare C code (DIGASO) and Matlab implementation of fbou
    % Flag implemented so we don't have to uncomment it every time
    if isfield(app,'check_fbou_flg')
        if (i == app.check_fbou_flg) && (size(fh,1)>0)      % Set app.check_fbou_flg to the index of the boundary you want to check

            nch=ncu;
            ncd=3;
            [ng,nc] = size(udgg);
            arg_dig = cell2mat(arg);
            ngb = ngf*nfb;

            fh_dig = zeros(size(fh));
            fh_udg_dig = zeros(size(dfh_dudg));
            fh_uh_dig = zeros(size(dfh_duh));

            % Flatten all arrays
            nl_dig = reshape(nlgb,1,[]);
            pg_dig = reshape(pgb,1,[]);
            udgg_dig = reshape(udgb,1,[]);
            uhg_dig = reshape(uhgb,1,[]);
            fh_dig = reshape(fh_dig,1,[]);
            fh_udg_dig = reshape(fh_udg_dig,1,[]);
            fh_uh_dig = reshape(fh_uh_dig,1,[]);
            bv_dig = reshape(bv,1,[]);

            a1=zeros(1,size(fh_dig,2));
            a2=zeros(1,size(fh_udg_dig,2));
            a3=zeros(1,size(fh_uh_dig,2));
            a4=zeros(1,size(pg_dig,2));
            a5=zeros(1,size(udgg_dig,2));
            a6=zeros(1,size(uhg_dig,2));
            a7=zeros(1,size(nl_dig,2));
            a8=zeros(1,size(bv_dig,2));
            a9=zeros(1,size(arg_dig,2));

            str = sprintf('Generating MEX file for fbou, boundary %d ...', i);
            disp(str)
            codegen fbou_digaso_matlab -args {a1,a2,a3,a4,a5,a6,a7,a8,a9,1,1,1,1,1,1,1} ../digaso_compare/fbou.c

            [fh_dig, fh_udg_dig, fh_uh_dig] = fbou_digaso_matlab_mex(fh_dig, fh_udg_dig, fh_uh_dig, pg_dig, udgg_dig, uhg_dig, nl_dig, bv_dig, arg_dig, time, ib, ngb, nc, nch, nd, ncd);

            % Reshape all arrays
            fh_dig = reshape(fh_dig, size(fh));
            fh_udg_dig = reshape(fh_udg_dig, size(dfh_dudg));
            fh_uh_dig = reshape(fh_uh_dig, size(dfh_duh));

            str=sprintf('Linf norm in fh: %d', max(max(fh_dig-fh)));
            disp(str)
            str=sprintf('Linf norm in fh_udg: %d', max(max(max(fh_udg_dig-dfh_dudg))));
            disp(str)
            str=sprintf('Linf norm in fh_uh: %d', max(max(max(fh_uh_dig-dfh_duh))));
            disp(str)

            error('Done checking fbou...')
        end
    end

    BH(im,:) = fh;
    BH_udg(im,:,:) = dfh_dudg;
    BH_uh(im,:,:) = dfh_duh;              
    
    if adjoint==1
        ia = bcd(i);        
        bu = repmat(bcv(i,:),[ngf*nfb 1]);            
        [~,djh_dudg,djh_duh] = fbou_adjoint(ia,bu,nlgb,pgb,udgb,uhgb,arg,time);                                               
        JUDG(im,:) = djh_dudg;
        JH(im,:) = djh_duh;
    end
end

if adjoint==1
    wrk = reshape(bsxfun(@times,JH,jac(bn)),[ngf nbf*nch]);
    Jh = zeros(npf*nfe*ne,nch);
    Jh(an,:) = reshape(shapfg*wrk,[npf*nbf nch]);
    Jh = permute(reshape(Jh,[npf*nfe ne nch]), [3 1 2]);

    wrk = reshape(bsxfun(@times,JUDG,jac(bn)),[ngf nbf*nc]);
    Jun = zeros(npf*nfe*ne,nc);
    Jun(an,:) = reshape(shapfg*wrk,[npf*nbf nc]);
    Jun = reshape(Jun,[npf nfe ne nc]);
    Juq = zeros(npv, 1, ne, nc);
    for is = 1:nfe
       IP = perm(:,is);
       Juq(IP,1,:,:) = Juq(IP,1,:,:) + Jun(:,is,:,:); 
    end
    Juq = reshape(Juq,[npv ne nc]);
    Ju  = permute(Ju + Juq(:,:,1:ncu),[1 3 2]);
    Jq  = permute(Jq + Juq(:,:,ncu+1:end),[1 3 2]);
    %if nc>ncu
    %    Jq  = reshape(Jq, [npv ncu nd ne]);
    %end
else
    Ju = [];
    Jq = [];
    Jh = [];    
end

wrk = reshape(bsxfun(@times,FH,jac),[ngf nfe*ne*nch]);
Run = reshape(shapfg*wrk,[npf nfe ne nch]);

% jac
% reshape(FH,[ngf*nfe ne*nch])
% reshape(bsxfun(@times,FH,jac),[ngf*nfe ne*nch])
% reshape(Run,[npf*nfe ne*nch])
% pause

wrk = reshape(bsxfun(@times,FH_udg,jac),[ngf nfe*ne*nch*nc]);
BDn = reshape(shapfgdotshapfc*wrk,[npf npf nfe ne nch nc]);
wrb = wrk;

wrk = reshape(bsxfun(@times,FH_uh,jac),[ngf nfe*ne*nch*nch]);
Fn = reshape(shapfgdotshapfc*wrk,[npf npf nfe ne nch nch]);
   
Rut = zeros(npv,1,ne,ncu);
BDt = zeros(npv,npv,1,ne,ncu,nc);
F = zeros(npv,npf,nfe,ne,ncu,nch); 
for is=1:nfe  % Adding face contributions - vector dependencies avoided
    IP = perm(:,is);
    Rut(IP,1,:,:) = Rut(IP,1,:,:) + Run(:,is,:,:);
    BDt(IP,IP,1,:,:,:) = BDt(IP,IP,1,:,:,:) + BDn(:,:,is,:,:,:);
    F(IP,:,is,:,:,:) = F(IP,:,is,:,:,:) + Fn(:,:,is,:,:,:);
end
Ru = permute(Ru-reshape(Rut,[npv ne ncu]),[1 3 2]);
BD = permute(BD+reshape(BDt,[npv npv ne ncu nc]),[1 4 2 5 3]);
F = permute(reshape(F,[npv npf*nfe ne ncu nch]),[1 4 5 2 3]);

wrk = reshape(bsxfun(@times,BH,jac(bn)),[ngf nbf*nch]);
Run = reshape(Run,[npf*nfe*ne nch]);
Run(an,:) = reshape(shapfg*wrk,[npf*nbf nch]);

Rh  = -permute(reshape(Run,[npf*nfe ne nch]), [3 1 2]);


wrk = reshape(bsxfun(@times,BH_udg,jac(bn)),[ngf nbf*nch*nc]);
Gt = reshape(BDn,[npf npf*nfe*ne nch nc]);
Gt(:,an,:,:) = reshape(shapfgdotshapfc*wrk,[npf npf*nbf nch nc]);
Gt = reshape(Gt,[npf npf nfe ne nch nc]);

wrk = reshape(bsxfun(@times,BH_uh,jac(bn)),[ngf nbf*nch*nch]);
Ht = reshape(Fn,[npf npf*nfe*ne nch nch]);
Ht(:,an,:,:) = reshape(shapfgdotshapfc*wrk,[npf npf*nbf nch nch]);
Ht  = reshape(Ht,[npf npf nfe ne nch nch]);

G   = zeros(npv,npf,nfe,ne,nch,nc); 
H   = zeros(npf*nfe,npf*nfe,1,ne,nch,nch);
for is=1:nfe  
    IP = (is-1)*npf+1:is*npf;
    G(perm(:,is),:,is,:,:,:) = G(perm(:,is),:,is,:,:,:) + Gt(:,:,is,:,:,:);
    H(IP,IP,1,:,:,:) = H(IP,IP,1,:,:,:) + Ht(:,:,is,:,:,:);
end
G  = permute(reshape(G,[npv npf*nfe ne nch nc]),[4 2 1 5 3]);
H  = permute(reshape(H,[npf*nfe npf*nfe ne nch nch]),[4 1 5 2 3]);
% t  = reshape(H(:,:,:,:,1),[nch*npf*nfe nch*npf*nfe])
% pause



% Read in arrays from digaso
% disp('DIGASO')
% % Ru, Rq, Rh, BD, F, E, G, H, Ju, Jq, Jh, wrb
% Rq_dig = readbin('/Users/saustin/Documents/digaso/problem/streamer/discharge_model/Rq.bin');
% Ru_dig = readbin('/Users/saustin/Documents/digaso/problem/streamer/discharge_model/Ru.bin');
% Rh_dig = readbin('/Users/saustin/Documents/digaso/problem/streamer/discharge_model/Rh.bin');
% BD_dig = readbin('/Users/saustin/Documents/digaso/problem/streamer/discharge_model/BD_mat.bin');
% F_dig = 
% E_dig = 
% G_dig = 
% H_dig = 
% Ju_dig = 
% Jq_dig = 
% Jh_dig = 
% wrb_dig = 
% 
% BD_dig = permute(BD_dig,[1,3,4,2]);
% % npv npv ne ncu nc
% 
% M_matlab = M(:,:,1);
% disp('MATLAB')
% M_matlab = chol(M_matlab);
% C_matlab = C(:,:,:,1);
% Rq_mat = Rq(:,:,:,1);
% Ru_mat = Ru(:,1,:);
% BD_mat = BD(:,:,1,:,:);
% 
% disp('M')
% max(abs(M_dig(:) - M_matlab(:)))
% disp('C')
% max(abs(C_dig(:) - C_matlab(:)))
% disp('Rq')
% max(abs(Rq_dig(:) - Rq_mat(:)))
% disp('Ru')
% max(abs(Ru_dig(:) - Ru_mat(:)))
% disp('BD')
% max(abs(BD_dig(:) - BD_mat(:)))
% disp('dgnodes')
% dgnodes = dgnodes(:,1,:);
% max(abs(dgnodes_dig(:) - dgnodes(:)))
% error('Done...')