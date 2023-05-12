function udg = ucg2udg(ucg, cgelcon)

[npe,ne] = size(cgelcon);
ncu = size(ucg,2);

udg = zeros(npe,ncu,ne);
for i = 1:ne
    udg(:,:,i) = ucg(cgelcon(:,i),:);
end

