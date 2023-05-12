function Ru = cuda_Ru_blockelem(Ru, xdgg, udgg, sg, Xx, jac, shapvg, param, time, fc_u, tdep, ne, ng, np, ncu, nd, nc)
% GPU implementation of volume integrals 

sg = reshape(sg,[ng*ne ncu nd]);
Xx = reshape(Xx,[ng*ne nd nd]);

% source at Gauss points
sg(:,:,1) = cuda_source(sg(:,:,1), xdgg, udgg, param, time); 
% Update source term for time-dependent problems
if tdep                
    sg(:,:,1) = sg(:,:,1) - udgg(:,1:ncu)*fc_u;    
end
for i = 1:ncu
    sg(:,i,1) = sg(:,i,1).*jac(:);
end

Ru(:,1:ncu) = cuda_gauss2node(Ru(:,1:ncu), sg(:,1:ncu,1), shapvg(:,:,1), ng, np, ne*ncu);

f = cuda_flux( xdgg, udgg, param, time);
f = reshape(f,[ng*ne ncu nd]);
for i=1:nd
    for m=1:ncu
        sg(:,m,i) = f(:,m,1).*Xx(:,1,i);
    end
    for j=2:nd
        for m=1:ncu
            sg(:,m,i) = sg(:,m,i) + f(:,m,j).*Xx(:,j,i);
        end
    end    
end

for i=1:nd
    Ru(:,1:ncu) = Ru(:,1:ncu) + cuda_gauss2node(Ru(:,1:ncu), sg(:,:,i), shapvg(:,:,i+1), ng, np, ne*ncu);
end

% storage: jac, xdgg(ngv,ne,nq), Xx(ngv,ne,nd,nd), udgg(ngv,ne,nc), sg, f, Ru
% Jg share memory with f(ngv,ne,max(ncu,nd),nd)
% pn, un and sn share memory with Ru(npv,ne,max(ncu,nq))

