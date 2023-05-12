function un = cuda_gauss2node(un, ug, shapg, ng, np, nn)

ncu = size(un,2);
un(:,1:ncu) = reshape(reshape(shapg,[np ng])*reshape(ug,[ng nn]),[],ncu);

