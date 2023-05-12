function [c, u] = uaverage(Wa, Ua, Pa, elcon)

c = sum(Wa(:).*Pa(:));

[~,nt] = size(Wa);
n = max(elcon(:));
u = zeros(1,n);
for i = 1:nt
    imap = elcon(:,i);  
    tm = Pa(:,i)'*Ua(:,:,i);
    u(imap) = u(imap) + tm;
end

% 1 * npe x npe * (npf*nfe) -> npf*nfe*ne -> nf*npf 

