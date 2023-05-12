function Rq = cuda_Rq_intface(Rq, Rqf, fhg, uhg, nlg1, jac1, shapfg, f2e, nf, ngf, npf, ncu)
% GPU implementation of interior face integrals 

for j = 1:nd
    for i = 1:ncu     
        fhg(:,i,j) = uhg(:,i).*(nlg1(:,j).*jac1(:));        
    end    
end

Rqf = cuda_gauss2node(Rqf, fhg, shapfg, ngf, npf, nf*ncu*nd, 1);

for i=1:npf*nf
    k1 = f2e(1,i);
    k2 = f2e(2,i);
    for m = 1:nd
        for j=1:ncu
            Rq(k1,j,m) = Rq(k1,j,m) - Rqf(i,j,m);  
            Rq(k2,j,m) = Rq(k2,j,m) + Rqf(i,j,m);  
        end
    end
end
    
end


