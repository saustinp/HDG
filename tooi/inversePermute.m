function B = inversePermute(A,order,maxDim)
l = 1:maxDim;

l(order) = 0;
r = sort(l(l>0));

newOrder = [order r];

if numel(newOrder) == 1
    newOrder = [1 2];
end

B = ipermute(A,newOrder);



% maxIndex = max(order);
% 
% inverseorder(order) = 1:maxIndex;   % Inverse permutation order
% a = permute(b,inverseorder);
% 
% i1 i2 i3 i4 i5
% 5  3
% 
% 1 --> 5
% 2 --> 3
% 3 --> 1
% 4 --> 2
% 5 --> 4
% 6 --> 6
% 
% l(order) = 0
% 
% sor(l(l > 0))
% 
% [order 
% 
% i1 1
% i2 2
% i3 3
% i4 4
% i5 5



