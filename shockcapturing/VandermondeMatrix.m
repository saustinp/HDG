function A = VandermondeMatrix(porder, plocal, elemtype)

if elemtype==0 % simplex elements
    A=koornwinder(plocal,porder);               % Vandermonde matrix
elseif elemtype==1           % tensor product elements
    A=tensorproduct(plocal,porder);             % Vandermonde matrix
end
   
