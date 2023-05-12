function M = mod(X,Y)
Y = mArray(Y);
M = mArray();
M.theArray = mod(X.theArray,Y.theArray);
end

