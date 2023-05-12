function u = getu(Wa, Ua, uhat, elcon)

[npv,nt] = size(Wa);
u = zeros(npv,1,nt);
for i = 1:nt
    imap = elcon(:,i);    
    u(:,1,i) = Wa(:,i) + Ua(:,:,i)*uhat(imap);
end

