function B = mPermute(A,indices)
maxIndices = max(ndims(A),max(indices));

l = 1:maxIndices;
l(indices) = 0;
r = sort(l(l>0));

newIndices = [indices r];

B = permute(A,newIndices);