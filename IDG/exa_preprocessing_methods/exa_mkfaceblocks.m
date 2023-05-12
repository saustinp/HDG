function [nm,nb] = exa_mkfaceblocks(mf,bcm,ns)

    if nargin<3
        ns = 4096; % default number of faces per block
    end
    
    nm = [];
    for i = 1:length(mf)-1
        nf = mf(i+1)-mf(i);    
        nmf = exa_mkelemblocks(nf,ns);
        tm = mf(i)+nmf;
        tm(3,:) = bcm(i);      
        nm = [nm tm];    
    end
    nb = size(nm,2);