function Ru = cuda_Ru_face_grad(Ru, udg, xdg, uinf, tempn, tempg, shapfn, shapfg, shapfgdotshapfc, param, time, f2e, nn, nmf,mesh)
    %Ruq = cuda_Ru_face(Ruq, permute(UDG,[1 3 2]), permute(mesh.dgnodes,[1 3 2]), uinf(1), tempn, tempg, shapfn, shapfg, param, time, f2e, nn, nmf);
    
    nd  = nn(1);
    nc  = nn(2);
    ncu = nn(3);
    ncx = nn(4);
    ne  = nn(5); % total number of elemens
    % nbe = nn(6); % number of element blocks
    % %blksze = nn(7); % maximum number of elements per blocks
    npe = nn(8);
    % nge = nn(9);
    %nf  = nn(10); % total number of faces
    nbf = nn(11); % number of face blocks
    %blkszf = nn(12); % maximum number of faces per blocks
    npf = nn(13);
    ngf = nn(14);
    Ru2 = Ru;
    
    for j=1:nbf(1)
        f1 = nmf(1,j);
        f2 = nmf(2,j);
        ib = nmf(3,j);
        ns = f2-f1+1;
        
        n1 = npf*ns;
        n2 = ngf*ns;
        tempn(1:n1,1:ncx) = cuda_getfacenodes(tempn(1:n1,1:ncx),xdg,f2e,npf,ncx,npe,f1,f2,1);
                
        % ncx+nd+1+nd*(nd-1)
        ind1 = 1:ncx; % xhg    
        ind2 = (ncx+1):(ncx+nd); % nlg
        ind3 = (ncx+nd+1); % jacg
        ind4 = (ncx+nd+2):(ncx+nd+1+nd*(nd-1)); % Jg
        [tempg(1:n2,ind1), tempg(1:n2,ind2), tempg(1:n2,ind3)] = cuda_facegeom(tempg(1:n2,ind1),... 
           tempg(1:n2,ind2), tempg(1:n2,ind3), tempg(1:n2,ind4), tempn(1:n1,1:ncx), shapfn, ns, ngf, npf, nd, ncx);
       
        tempn(1:n1,1:ncu) = cuda_getfacenodes(tempn(1:n1,1:ncu),udg,f2e,npf,ncu,npe,f1,f2,0);        
        ind5 = (ncx+nd+2):(ncx+nd+1+ncu); %uhg
        tempg(1:n2,ind5) = cuda_node2gauss(tempg(1:n2,ind5), tempn(1:n1,1:ncu), shapfn, ngf, npf, ns*ncu);            
    
        tempn(1:n1,1:nc) = cuda_getfacenodes(tempn(1:n1,1:nc),udg,f2e,npf,nc,npe,f1,f2,1); 
        ind6 = (ncx+nd+2+ncu):(ncx+nd+1+ncu+nc); %udgg1
        tempg(1:n2,ind6) = cuda_node2gauss(tempg(1:n2,ind6), tempn(1:n1,1:nc), shapfn, ngf, npf, ns*nc);            
        
        tempn(1:n1,1:nc) = cuda_getfacenodes(tempn(1:n1,1:nc),udg,f2e,npf,nc,npe,f1,f2,2); 
        ind7 = (ncx+nd+2+ncu+nc):(ncx+nd+1+ncu+nc+nc); %udgg2
        tempg(1:n2,ind7) = cuda_node2gauss(tempg(1:n2,ind7), tempn(1:n1,1:nc), shapfn, ngf, npf, ns*nc);            
          
        % ncx+nd+1+max(2*nc+2*ncu,nd*(nd-1))            
        ind8 = (ncx+nd+2+ncu+nc+nc):(ncx+nd+1+ncu+nc+nc+ncu);
        tempn(1:n1,1:ncu) = cuda_Ru_blockface_grad(tempn(1:n1,1:ncu), tempg(1:n2,ind8), tempg(1:n2,ind6), tempg(1:n2,ind7), ...
             tempg(1:n2,ind5), uinf, tempg(1:n2,ind1), tempg(1:n2,ind2), tempg(1:n2,ind3), shapfg, shapfgdotshapfc, ... 
             param, time, ib, ns, ngf, npf, ncu, nd, nc);   
    
        eps = 1e-6;
        [Re, Re_dup, Re_dum] = cuda_Ru_blockface_grad(tempn(1:n1,1:ncu), tempg(1:n2,ind8), tempg(1:n2,ind6), tempg(1:n2,ind7), ...
             tempg(1:n2,ind5), uinf, tempg(1:n2,ind1), tempg(1:n2,ind2), tempg(1:n2,ind3), shapfg,shapfgdotshapfc, ... 
             param, time, ib, ns, ngf, npf, ncu, nd, nc);   
    
    %TODO: finite difference check here...
    
    %     Re_up_eps = cuda_Ru_blockface_grad(tempn(1:n1,1:ncu), tempg(1:n2,ind8), tempg(1:n2,ind6) * (1 + eps), tempg(1:n2,ind7), ...
    %          tempg(1:n2,ind5), uinf, tempg(1:n2,ind1), tempg(1:n2,ind2), tempg(1:n2,ind3), shapfg, shapfgdotshapfc,... 
    %          param, time, ib, ns, ngf, npf, ncu, nd, nc); 
    
    %     Re_um_eps = cuda_Ru_blockface_grad(tempn(1:n1,1:ncu), tempg(1:n2,ind8), tempg(1:n2,ind6), tempg(1:n2,ind7)*(1 + eps), ...
    %          tempg(1:n2,ind5), uinf, tempg(1:n2,ind1), tempg(1:n2,ind2), tempg(1:n2,ind3), shapfg, shapfgdotshapfc,... 
    %          param, time, ib, ns, ngf, npf, ncu, nd, nc); 
    
    %     grad_up = (Re_up_eps - Re) / eps;
    %     grad_um = (Re_um_eps - Re) / eps;
    %      
%         JRu = zeros([size(Ru), size(Ru)]); %Should probably be sparsified...
        Ru = cuda_putfacenodes(Ru,tempn(1:n1,1:ncu),f2e,npf,ncu,npe,f1,f2,ib); 
%         JRu = cuda_putfacenodes_Jac(JRu,Re_dup,Re_dum,f2e,npf,ncu,npe,f1,f2,ib); 
    end
    Ru2 = cuda_putfacenodes_nodal_based(Ru2, tempn(1:n1,1:ncu), mesh.rowe2f1, mesh.cole2f1, mesh.ent2ind1, mesh.rowe2f2, mesh.cole2f2, mesh.ent2ind2, npf, npe, nc, 1, ne, 0)

    
    