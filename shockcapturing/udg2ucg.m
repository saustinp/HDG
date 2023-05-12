function vcg = udg2ucg(udg, cgent2dgent, rowent2elem)

[npe,ncu,ne] = size(udg);
udg = reshape(permute(udg, [1 3 2]),[npe*ne ncu]);

nent = length(rowent2elem)-1;
vcg = zeros(nent, ncu);
for i = 1:nent
    dgn = cgent2dgent((rowent2elem(i)+1):rowent2elem(i+1));
    vcg(i,:) = mean(udg(dgn,:),1);
end

