function uh = cuda_getfacenodes(uh,udg,f2e,npf,ncu,npe,f1,f2,opts)

nf = f2-f1+1;
ndf = npf*nf;
if opts==0
    for i=1:ndf
        m = npf*(f1-1)+i;
        k1 = f2e(1,m); m1 = rem(k1-1,npe)+1; n1 = (k1-m1)/npe+1;
        k2 = f2e(2,m); m2 = rem(k2-1,npe)+1; n2 = (k2-m2)/npe+1;   
        for j=1:ncu
            uh(i+(j-1)*ndf) = 0.5*(udg(m1+(j-1)*npe+(n1-1)*npe*ncu)+udg(m2+(j-1)*npe+(n2-1)*npe*ncu));                
        end
    end
elseif opts==1    
    for i=1:ndf        
        m = npf*(f1-1)+i;
        k1 = f2e(1,m); m1 = rem(k1-1,npe)+1; n1 = (k1-m1)/npe+1;        
        for j=1:ncu
            uh(i+(j-1)*ndf) = udg(m1+(j-1)*npe+(n1-1)*npe*ncu);
        end
    end    
elseif opts==2
    for i=1:ndf
        m = npf*(f1-1)+i;
        k2 = f2e(2,m); m2 = rem(k2-1,npe)+1; n2 = (k2-m2)/npe+1;       
        for j=1:ncu
            uh(i+(j-1)*ndf) = udg(m2+(j-1)*npe+(n2-1)*npe*ncu);
        end
    end        
end

% function uh = cuda_getfacenodes(uh,udg,f2e,npf,ncu,nde,f1,f2,opts)
% 
% nf = f2-f1+1;
% ndf = npf*nf;
% if opts==0
%     for i=1:ndf
%         m = npf*(f1-1)+i;
%         k1 = f2e(1,m); % m1 = rem(k1-1,npe)+1 and n1 = (k1-m1)/npe+1;
%         k2 = f2e(2,m);    
%         for j=1:ncu
%             uh(i+(j-1)*ndf) = 0.5*(udg(k1+(j-1)*nde)+udg(k2+(j-1)*nde));                
%         end
%     end
% elseif opts==1    
%     for i=1:ndf        
%         m = npf*(f1-1)+i;
%         k1 = f2e(1,m);        
%         for j=1:ncu
%             uh(i+(j-1)*ndf) = udg(k1+(j-1)*nde); 
%         end
%     end    
% elseif opts==2
%     for i=1:ndf
%         m = npf*(f1-1)+i;
%         k2 = f2e(2,m);    
%         for j=1:ncu
%             uh(i+(j-1)*ndf) = udg(k2+(j-1)*nde);                
%         end
%     end        
% end
