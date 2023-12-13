function q = getq(MiB, MiC, u, uhat, elcon)

[npv,~,nd,nt] = size(MiC);

q = zeros(npv,nd,nt);
for i = 1:nt
    imap = elcon(:,i);    
    for l = 1:nd      
      q(:,l,i) = MiC(:,:,l,i)*uhat(imap) - MiB(:,:,l,i)*u(:,1,i);     
    end    
end

