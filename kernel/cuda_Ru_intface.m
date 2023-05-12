function Ru = cuda_Ru_intface(Ru, Ruf, fhg, udg1, udg2, uhg, pg1, nlg1, jac1, shapfg, param, time, f2e, nf, ngf, npf, ncu)
% GPU implementation of interior face integrals 

fhg = fhg + dgfhat(nlg1, pg1, udg1, udg2, uhg, param, time);
for i = 1:ncu
    fhg(:,i) = fhg(:,i).*jac1(:);
end

Ruf = cuda_gauss2node(Ruf, fhg, shapfg, ngf, npf, nf*ncu, 1);

for i=1:npf*nf
    k1 = f2e(i,1);
    k2 = f2e(i,2);
    for j=1:ncu
        Ru(k1,j) = Ru(k1,j) - Ruf(i,j);  
        Ru(k2,j) = Ru(k2,j) + Ruf(i,j);  
    end
end

