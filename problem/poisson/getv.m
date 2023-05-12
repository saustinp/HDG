function v = getv(MiD, MiE, u, uhat, elcon)

[npv,~,nd,~,nt] = size(MiD);

v = zeros(npv,nd,nd,nt);
for i = 1:nt
    imap = elcon(:,i);    
    for l = 1:nd
      for m = 1:nd
        v(:,l,m,i) = MiD(:,:,l,m,i)*u(:,1,i) + MiE(:,:,l,m,i)*uhat(imap);     
      end
    end    
end

