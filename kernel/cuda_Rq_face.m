function [Rq,xhg,nlg,jcg,ug1,ug2,uhg] = cuda_Rq_face(Rq, udg, xdg, uh, uinf, tempn, tempg, shapfn, shapfg, param, time, f2e, nn, nmf)
%function Ru = cuda_Ru_face(Ru, udg, xdg, uinf, tempn, tempg, shapfn, shapfg, param, time, f2e, nn, nmf)

nd  = nn(1);
nc  = nn(2);
ncu = nn(3);
ncx = nn(4);
ne  = nn(5); % total number of elemens
% nbe = nn(6); % number of element blocks
% %blksze = nn(7); % maximum number of elements per blocks
npe = nn(8);
% nge = nn(9);
nf  = nn(10); % total number of faces
nbf = nn(11); % number of face blocks
%blkszf = nn(12); % maximum number of faces per blocks
npf = nn(13);
ngf = nn(14);

for j=1:size(nmf,2)
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
    xhg = tempg(1:n2,ind1);
    nlg = tempg(1:n2,ind2);
    jcg = tempg(1:n2,ind3);  
%     tempn(1:n1,1:ncu) = cuda_getfacenodes(tempn(1:n1,1:ncu),udg,f2e,npf,ncu,ne*npe,f1,f2,0);        
    %ind5 = (ncx+nd+2):(ncx+nd+1+ncu); %uhg
%     tempg(1:n2,ind5) = cuda_node2gauss(tempg(1:n2,ind5), tempn(1:n1,1:ncu), shapfn, ngf, npf, ns*ncu);            

    tempn(1:n1,1:nc) = cuda_getfacenodes(tempn(1:n1,1:nc),udg,f2e,npf,nc,npe,f1,f2,1); 
    ind6 = (ncx+nd+2):(ncx+nd+1+nc); %udgg1
    tempg(1:n2,ind6) = cuda_node2gauss(tempg(1:n2,ind6), tempn(1:n1,1:nc), shapfn, ngf, npf, ns*nc);            
    ug1 = tempg(1:n2,ind6);
    
    tempn(1:n1,1:nc) = cuda_getfacenodes(tempn(1:n1,1:nc),udg,f2e,npf,nc,npe,f1,f2,2); 
    ind7 = (ncx+nd+2+nc):(ncx+nd+1+nc+nc); %udgg2
    tempg(1:n2,ind7) = cuda_node2gauss(tempg(1:n2,ind7), tempn(1:n1,1:nc), shapfn, ngf, npf, ns*nc);            
    ug2 = tempg(1:n2,ind7);
     
    % ncx+nd+1+2*nc+ncu*nd            
    ind8 = (ncx+nd+2+nc+nc):(ncx+nd+1+nc+nc+ncu*nd); %uhg
    tempn(1:n1,1:ncu) = cuda_getelemnodes(tempn(1:n1,1:ncu), uh, npf, nf, ncu, f1, f2);    
    tempg(1:n2,ind8) = cuda_node2gauss(tempg(1:n2,ind8), tempn(1:n1,1:ncu), shapfn, ngf, npf, ns*ncu);                
    uhg = tempg(1:n2,ind8);
    
    tempn(1:n1,1:ncu*nd) = cuda_Rq_blockface(tempn(1:n1,1:ncu*nd), tempg(1:n2,ind8), tempg(1:n2,ind6), tempg(1:n2,ind7), ...
         uinf, tempg(1:n2,ind1), tempg(1:n2,ind2), tempg(1:n2,ind3), shapfg, ... 
         param, time, ib, ns, ngf, npf, ncu, nd);        
    
    Rq = cuda_putfacenodes(Rq,tempn(1:n1,1:ncu*nd),f2e,npf,ncu*nd,npe,f1,f2,ib);
end


