function Rq = cuda_Rq_bouface(Rq, Rqf, fhg, uhg, nlg1, jac1, shapfg, f2e, ib, nf, ngf, npf, ncu)
% GPU implementation of boundary face integrals 

uhg = dgubou(ib, uhg, uinf, nlg1, pg1, udg1, f2e);
for j = 1:nd
    for i = 1:ncu     
        fhg(:,i,j) = uhg(:,i).*(nlg1(:,j).*jac1(:));        
    end    
end

Rqf = cuda_gauss2node(Rqf, fhg, shapfg, ngf, npf, nf*ncu*nd, 1);

for i=1:npf*nf
    k1 = f2e(i,1);    
    for m = 1:nd
        for j=1:ncu
            Rq(k1,j,m) = Rq(k1,j,m) - Rqf(i,j,m);              
        end
    end
end
    
end


