function [Rq,xdgg,Xxg,jacg,udgg] = cuda_Rq_elem(Rq, udg, xdg, tempn, tempg, shapen, shapeg, nn, nm)

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
    xdgg = tempg(1:n2,ind1);
    Xxg = tempg(1:n2,ind2);
    jacg = tempg(1:n2,ind3);    
    
    tempn(1:n1,1:nc) = cuda_getelemnodes(tempn(1:n1,1:nc), udg, npe, ne, nc, e1, e2);    
    ind5 = (ncx+nd*nd+2):(ncx+nd*nd+1+nc); % udgg
    tempg(1:n2,ind5) = cuda_node2gauss(tempg(1:n2,ind5), tempn(1:n1,1:nc), shapen, nge, npe, ns*nc);
    udgg = tempg(1:n2,ind5);
    
    ind4 = (ncx+nd*nd+1+nc+1):(ncx+nd*nd+1+nc+ncu*nd); % 
    tempn(1:n1,1:ncu*nd) = cuda_Rq_blockelem(tempn(1:n1,1:ncu*nd), tempg(1:n2,ind5), tempg(1:n2,ind4), tempg(1:n2,ind2), ...
        shapeg, ns, nge, npe, ncu, nd);
    
    %Rq = tempn(1:n1,1:ncu*nd);
    Rq = cuda_putelemnodes(Rq, tempn(1:n1,1:ncu*nd), npe, ne, ncu*nd, e1, e2);
end

