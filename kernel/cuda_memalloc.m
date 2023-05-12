function [Ruq,tempn,tempg] = cuda_memalloc(nn)

nd  = nn(1);
nc  = nn(2);
ncu = nn(3);
ncx = nn(4);
ne  = nn(5); % total number of elemens
%nbe = nn(6); % number of element blocks
blksze = nn(7); % maximum number of elements per blocks
npe = nn(8);
nge = nn(9);

n1 = max(nc+ncx+1,ncu*nd)+nd*nd+ncu*nd;
n2 = ncx+1+nd*nd+max(nd*nd,nc+ncu*nd);
tempg = zeros(nge*blksze,max([n1,n2]));

tempn = zeros(npe*blksze,max([ncx,nc]));

Ruq = zeros(npe,ncu*nd,ne);


% blksze = 512;
% [nm,nb] = cuda_mkblocks(ne,blksze);
% 
% ng = nge*blksze;
% tempg = cuda_memalloc(ng,nc,ncu,ncx,nd);
% 
% np = npe*blksze;
% tempn = zeros(np,max(ncx,nc));
% 
% % initialize
% Rq = zeros(npe, ne, nc);


% udgg = zeros(nge*blksze, nc);
% xdgg = zeros(nge*blksze, ncx);
% jacg = zeros(nge*blksze, 1);    
% Xxg = zeros(nge*blksze, nd, nd);    
% sg = zeros(nge*blksze, ncu, nd);
% fg = zeros(nge*blksze, ncu, nd);
