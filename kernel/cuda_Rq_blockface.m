function Rqf = cuda_Rq_blockface(Rqf, uhg, udg1, udg2, uinf, xdg1, nlg1, jac1, shapfg, param, time, ib, nf, ngf, npf, ncu, nd)

% uhg = cuda_uhat(ib, udg1, udg2, uhg, uinf, xdg1, nlg1, jac1, param, time, ngf*nf, ncu, nd);
% fhg = reshape(fhg,[ngf*nf,ncu,nd]);
% for j = 1:nd
%     for i = 1:ncu     
%         fhg(:,i,j) = uhg(:,i).*(nlg1(:,j).*jac1(:));        
%     end    
% end
% fhg = reshape(fhg,[ngf*nf,ncu*nd]);
uhg = reshape(uhg,[ngf*nf,ncu,nd]);
for j = nd:-1:1    
    for i = 1:ncu         
        uhg(:,i,j) = uhg(:,i,1).*nlg1(:,j).*jac1(:);
    end    
end    
uhg = reshape(uhg,[ngf*nf,ncu*nd]);

Rqf(:,1:ncu*nd) = cuda_gauss2node(Rqf(:,1:ncu*nd), uhg, shapfg(:,:,1), ngf, npf, nf*ncu*nd);
