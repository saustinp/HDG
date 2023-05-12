function fhg = cuda_fhat(ib, udg1, udg2, uhg, uinf, xdg1, nlg1, jac1, param, time, ng, ncu, nd, nc)
ib
if ib==0
    fhg = ldgfhat(nlg1, xdg1, udg1, udg2, uhg, param, time);
elseif ib>0
    fhg = ldgfbou(ib, uinf(1:ncu), nlg1, xdg1, udg1, uhg, param, time);
end

for i = 1:ncu
    fhg(:,i) = fhg(:,i).*jac1(:);
end

