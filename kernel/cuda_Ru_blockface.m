function Ruf = cuda_Ru_blockface(Ruf, fhg, udg1, udg2, uhg, uinf, xdg1, nlg1, jac1, shapfg, param, time, ib, nf, ngf, npf, ncu, nd, nc)
% GPU implementation of face integrals 

fhg(:,1:ncu) = cuda_fhat(ib, udg1, udg2, uhg, uinf, xdg1, nlg1, jac1, param, time, ngf*nf, ncu, nd, nc);

Ruf(:,1:ncu) = cuda_gauss2node(Ruf(:,1:ncu), fhg(:,1:ncu), shapfg(:,:,1), ngf, npf, nf*ncu);

%Rqf(:,1:ncu*nd) = cuda_gauss2node(Rqf(:,1:ncu*nd), uhg, shapfg(:,:,1), ngf, npf, nf*ncu*nd);
