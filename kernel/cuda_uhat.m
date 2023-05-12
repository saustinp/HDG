function uhg = cuda_uhat(ib, udg1, udg2, uhg, uinf, xdg1, nlg1, jac1, param, time, ng, ncu, nd, nc)

uhg = reshape(uhg,[ng,ncu,nd]);
if ib==0
    for i=1:ncu
        uhg(:,i,1) = (udg1(:,i) + udg2(:,i))/2 + sign(nlg1(:,1)+nlg1(:,2)).*(udg1(:,i) - udg2(:,i))/2;
    end
elseif ib>0
    uhg(:,1:ncu,1) = ldgubou(ib,uinf,nlg1,xdg1,udg1,param,time);
end
for j = nd:-1:1    
    for i = 1:ncu         
        uhg(:,i,j) = uhg(:,i,1).*nlg1(:,j).*jac1(:);
    end    
end    
uhg = reshape(uhg,[ng,ncu*nd]);


% if ib>0 
%     % boundary conditions
%     udg2 = ldgubou(ib,uinf,nlg1,xdg1,udg1,param,time);
% end
% 
% uhg = reshape(uhg,[ng,ncu,nd]);
% for j = 1:nd
%     for i = 1:ncu     
%         uhg(:,i,j) = 0.5*((udg1(:,i) + udg2(:,i)).*nlg1(:,j)+(udg1(:,i) - udg2(:,i))).*jac1(:);        
%     end    
% end    
% uhg = reshape(uhg,[ng,ncu*nd]);
