function ug = cuda_node2gauss(ug, un, shapt, ng, np, nn)

ncu = size(un,2);
ug(:,1:ncu) = reshape(reshape(shapt(:,:,1),[ng np])*reshape(un,[np nn]),[],ncu);
