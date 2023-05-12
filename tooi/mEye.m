function array = mEye(sizes,i,j)
    array = mArray.zeros(sizes,[i j]);

    nk = min(sizes(1),sizes(2));
    for k = 1:nk
        array(i,k,j,k) = 1;
    end            
end