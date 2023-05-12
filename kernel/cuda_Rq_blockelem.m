function Rq = cuda_Rq_blockelem(Rq, udgg, fqg, Xx, shapvg, ne, ng, np, ncu, nd)
% GPU implementation of volume integrals 

Xx = reshape(Xx,[ng*ne,nd,nd]);
Rq = reshape(Rq(:,1:ncu*nd),[np*ne,ncu,nd]);
for i=1:nd
    for m=1:ncu
        fqg(:,m,i) = udgg(:,m).*Xx(:,1,i);
    end
    for j=2:nd
        for m=1:ncu
            fqg(:,m,i) = fqg(:,m,i) + udgg(:,m).*Xx(:,j,i);
        end
    end    
    Rq(:,:,i) = cuda_gauss2node(Rq(:,:,i), fqg(:,1:ncu,i), shapvg(:,:,i+1), ng, np, ne*ncu);
end
Rq = reshape(Rq,[np*ne,ncu*nd]);





