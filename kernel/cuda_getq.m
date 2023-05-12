function q = cuda_getq(q, Mi, Rq, sdg, udg, xdg, uh, uinf, tempen, tempeg, tempfn, tempfg, shapen, shapeg, shapfn, shapfg, param, time, fc_q, f2e, nn, nme, nmf)

Rq = cuda_Rq_elem(Rq, udg, xdg, tempen, tempeg, shapen, shapeg, nn, nme);
Rq = cuda_Rq_face(Rq, udg, xdg, uh, uinf, tempfn, tempfg, shapfn, shapfg, param, time, f2e, nn, nmf);

nd = nn(1); 
ncu = nn(3);
ne = nn(5);
npe = nn(8);
ncq = ncu*nd;
Rq = reshape(Rq,[npe ncq ne]);
q = reshape(q,[npe ncq ne]);
for i = 1:ne    
    q(:,:,i) = Mi(:,:,i)*Rq(:,:,i);    
end
q = (1/fc_q)*(q + sdg(:,ncu+1:ncu+ncq,:));        

% permute Rq
% for i = 1:npe
%     for j = 1:ne        
%         for k = 1:ncq     
%             m = i + (k-1)*npe + (j-1)*npe*ne;
%             q(m) = Rq(i,j,k);            
%         end
%     end
% end
% for n = 0:npe*ne*ncq-1
%     i = rem(n,npe); %[0,1,...,npe-1]
%     j = rem(n,ncq); %[0,1,...,ncq-1]
%     k = rem(n,ne);  %[0,1,...,ne-1]
%     m = i + j*npe + k*npe*ne;
%     q(m+1) = Rq(n+1);
% end




