function [Ruf, Ruf_dudg1, Ruf_dudg2] = cuda_Ru_blockface_grad(Ruf, fhg, udg1, udg2, uhg, uinf, xdg1, nlg1, jac1, shapfg, shapfgdotshapfc, param, time, ib, nf, ngf, npf, ncu, nd, nc)
    % GPU implementation of face integrals 
    
    % fhg_udg1 = zeros(size)
    % 
    %TODO: Shoudl be preallocated
    fhg_udg1 = zeros(size(fhg,1),size(fhg,2),ncu);
    fhg_udg2 = zeros(size(fhg,1),size(fhg,2),ncu);
    
    [fhg(:,1:ncu), fhg_udg1(:,1:ncu,1:ncu), fhg_udg2(:,1:ncu,1:ncu)] = ...
        cuda_fhat_grad(ib, udg1, udg2, uhg, uinf, xdg1, nlg1, jac1, param, time, ngf*nf, ncu, nd, nc);
    
    Ruf(:,1:ncu) = cuda_gauss2node(Ruf(:,1:ncu), fhg(:,1:ncu), shapfg(:,:,1), ngf, npf, nf*ncu);
    
    Ruf_dudg1 = reshape(shapfgdotshapfc(:,:,1), [npf*npf ngf]) * reshape(fhg_udg1, [ngf ncu*ncu*nf]);
    Ruf_dudg2 = reshape(shapfgdotshapfc(:,:,1), [npf*npf ngf]) * reshape(fhg_udg2, [ngf ncu*ncu*nf]);
    
    Ruf_dudg1 = reshape(Ruf_dudg1, [npf ncu npf ncu nf]);
    Ruf_dudg2 = reshape(Ruf_dudg2, [npf ncu npf ncu nf]);
    
    
    %Rqf(:,1:ncu*nd) = cuda_gauss2node(Rqf(:,1:ncu*nd), uhg, shapfg(:,:,1), ngf, npf, nf*ncu*nd);
    