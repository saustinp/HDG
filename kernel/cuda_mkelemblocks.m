function [nm,nb] = cuda_mkelemblocks(ne,ns)

if nargin<2
    ns = 512; % default number of elements per block
end

if ne < ns
   ns = ne;
end

nb = ceil(ne/ns);          
nk = 1:ns:ne;
nm = [nk(1:end); [nk(2:end)-1,ne]];    
