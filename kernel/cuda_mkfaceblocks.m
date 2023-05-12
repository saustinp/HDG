function [nm,nb] = cuda_mkfaceblocks(nf,bcm,ns)

if nargin<3
    ns = 512; % default number of faces per block
end

nb = zeros(length(nf)-1,1);
nm = [];
for i = 1:length(nf)-1
    nfi = nf(i+1)-nf(i);    
        
    if nfi < ns
       n = nfi;
    else
       n = ns; 
    end

    nb(i) = ceil(nfi/n);          
    nk = 1:n:nfi;
    tm = nf(i)+[nk(1:end); [nk(2:end)-1,nfi]];        
    tm(3,:) = bcm(i);      
    nm = [nm tm];
end


