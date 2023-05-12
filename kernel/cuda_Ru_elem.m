function Ru = cuda_Ru_elem(Ru, udg, xdg, sdg, tempn, tempg, shapen, shapeg, param, time, fc_u, tdep, nn, nm)

nd  = nn(1);
nc  = nn(2);
ncu = nn(3);
ncx = nn(4);
ne  = nn(5); % total number of elemens
nbe = nn(6); % number of element blocks
%blksze = nn(7); % maximum number of elements per blocks
npe = nn(8);
nge = nn(9);

% compute volume integrals to form Rq
for j=1:nbe
    e1 = nm(1,j);
    e2 = nm(2,j);    
    ns = e2-e1+1;
    n1 = npe*ns;
    n2 = nge*ns;
                            
    tempn(1:n1,1:ncx) = cuda_getelemnodes(tempn(1:n1,1:ncx), xdg, npe, ne, ncx, e1, e2);    
    ind1 = 1:ncx; % xdgg    
    ind2 = (ncx+1):(ncx+nd*nd); % Xxg
    ind3 = (ncx+nd*nd+1); % jacg
    ind4 = (ncx+nd*nd+2):(ncx+nd*nd+1+nd*nd); % Jg
    [tempg(1:n2,ind1), tempg(1:n2,ind2), tempg(1:n2,ind3)] = cuda_elemgeom(tempg(1:n2,ind1), ...
        tempg(1:n2,ind2), tempg(1:n2,ind3), tempn(1:n1,1:ncx), tempg(1:n2,ind4), shapen, ns, nge, npe, nd, ncx);
    
    tempn(1:n1,1:nc) = cuda_getelemnodes(tempn(1:n1,1:nc), udg, npe, ne, nc, e1, e2);    
    ind4 = (ncx+nd*nd+2):(ncx+nd*nd+1+nc); % udgg
    tempg(1:n2,ind4) = cuda_node2gauss(tempg(1:n2,ind4), tempn(1:n1,1:nc), shapen, nge, npe, ns*nc);
        
    tempn(1:n1,1:ncu) = cuda_getelemnodes(tempn(1:n1,1:ncu), sdg, npe, ne, ncu, e1, e2);  
    ind5 = (ncx+nd*nd+1+nc+1):(ncx+nd*nd+1+nc+ncu*nd); % source term due to time derivative
    tempg(1:n2,ind5(1:ncu)) = cuda_node2gauss(tempg(1:n2,ind5(1:ncu)), tempn(1:n1,1:ncu), shapen, nge, npe, ns*ncu);
    
    tempn(1:n1,1:ncu) = cuda_Ru_blockelem(tempn(1:n1,1:ncu), tempg(1:n2,ind1), tempg(1:n2,ind4), tempg(1:n2,ind5), ...
        tempg(1:n2,ind2), tempg(1:n2,ind3), shapeg, param, time, fc_u, tdep, ns, nge, npe, ncu, nd);
    
    Ru = cuda_putelemnodes(Ru, tempn(1:n1,1:ncu), npe, ne, ncu, e1, e2);
end

