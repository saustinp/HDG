function Ru = cuda_Ru_bouface(Ru, Ruf, udg1, uinf, uhg, pg1, nlg1, jac1, shapfg, f2e, param, time, ib, nf, ngf, npf, ncu)
%function Ru = cuda_Ru_intface(Ru, Ruf, fhg, udg1, udg2, uhg, pg1, nlg1, jac1, shapfg, param, time, f2e, nf, ngf, npf, ncu)
% GPU implementation of interior face integrals 

fhg = dgfbou(ib, uinf, nlg1, pg1, udg1, uhg, param, time);
for i = 1:ncu
    fhg(:,i) = fhg(:,i).*jac1(:);
end

Ruf = cuda_gauss2node(Ruf, fhg, shapfg, ngf, npf, nf*ncu, 1);

for i=1:npf*nf
    k1 = f2e(i,1);
    for j=1:ncu
        Ru(k1,j) = Ru(k1,j) - Ruf(i,j);  
    end
end

