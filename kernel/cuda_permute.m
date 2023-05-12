function [A,dimA] = cuda_permute(B, dimB, order, d)
% CUDA implementation of matlab permutation: A = permute(reshape(B,dimB),order)
% Example:
% d = 4;
% order = randperm(d);
% dimB = randi(10,[1 d]);
% B = rand(dimB);
% [A,dimA] = cuda_permute(B, dimB, order, d);

dimA = dimB(order(1:d));
subA = zeros(1,d);
subB = zeros(1,d);
nB = prod(dimB(1:d));
A = zeros(nB,1);

lA = zeros(max(d-2,1),1);
lB = zeros(max(d-2,1),1);
lA(1) = dimA(1)*dimA(2);
lB(1) = dimB(1)*dimB(2);
for i = 2:d-2    
    lA(i) = lA(i-1)*dimA(i+1);    
    lB(i) = lB(i-1)*dimB(i+1);    
end

if d==2
    A = permute2d(A, B, dimA, dimB, subA, subB, order, nB);
elseif d==3
    A = permute3d(A, B, dimA, dimB, subA, subB, lA, lB, order, nB);    
elseif d==4
    A = permute4d(A, B, dimA, dimB, subA, subB, lA, lB, order, nB);        
elseif d==5
    A = permute5d(A, B, dimA, dimB, subA, subB, lA, lB, order, nB);            
elseif d==6
    A = permute6d(A, B, dimA, dimB, subA, subB, lA, lB, order, nB);                
elseif d==7
    A = permute7d(A, B, dimA, dimB, subA, subB, lA, lB, order, nB);         
else
    error('Array dimensions must be in the range [2, 7]');
end

function A = permute2d(A, B, dimA, dimB, subA, subB, order, nB)

% for subB1 = 1:dimB(1)
%     for subB2 = 1:dimB(2)        
%         subB(1) = subB1;
%         subB(2) = subB2; 
%         indB = subB(1) + (subB(2)-1)*dimB(1); 
%         subA(1) = subB(order(1));
%         subA(2) = subB(order(2));
%         indA = subA(1) + (subA(2)-1)*dimA(1);        
%         A(indA) = B(indB);
%     end
% end
for n = 0:nB-1
    subB(1) = rem(n,dimB(1))+1;
    subB(2) = (n+1-subB(1))/dimB(1)+1;     
    indB = subB(1) + (subB(2)-1)*dimB(1); 
    subA(1) = subB(order(1));
    subA(2) = subB(order(2));
    indA = subA(1) + (subA(2)-1)*dimA(1);        
    A(indA) = B(indB);
end


function A = permute3d(A, B, dimA, dimB, subA, subB, lA, lB, order, nB)

% for subB1 = 1:dimB(1)
%     for subB2 = 1:dimB(2)       
%         for subB3 = 1:dimB(3)       
%             subB(1) = subB1;
%             subB(2) = subB2; 
%             subB(3) = subB3; 
%             indB = subB(1) + (subB(2)-1)*dimB(1) + (subB(3)-1)*dimB(1)*dimB(2); 
%             subA(1) = subB(order(1));
%             subA(2) = subB(order(2));
%             subA(3) = subB(order(3));
%             indA = subA(1) + (subA(2)-1)*dimA(1) + (subA(3)-1)*dimA(1)*dimA(2);         
%             A(indA) = B(indB);
%         end
%     end
% end
for n = 0:nB-1
    subB(1) = rem(n,dimB(1))+1;
    k = rem(n,lB(1));        
    subB(2) = (k+1-subB(1))/dimB(1)+1;  
    subB(3) = (n+1-subB(1)-(subB(2)-1)*dimB(1))/(lB(1))+1;      
    indB = subB(1) + (subB(2)-1)*dimB(1) + (subB(3)-1)*lB(1); 
    subA(1) = subB(order(1));
    subA(2) = subB(order(2));
    subA(3) = subB(order(3));
    indA = subA(1) + (subA(2)-1)*dimA(1) + (subA(3)-1)*lA(1);             
    A(indA) = B(indB);
end

function A = permute4d(A, B, dimA, dimB, subA, subB, lA, lB, order, nB)

for n = 0:nB-1
    subB(1) = rem(n,dimB(1))+1;
    k = rem(n,lB(1));        
    subB(2) = (k+1-subB(1))/dimB(1)+1;  
    k = rem(n,lB(2));        
    subB(3) = (k+1-subB(1)-(subB(2)-1)*dimB(1))/(lB(1))+1;   
    subB(4) = (n+1-subB(1)-(subB(2)-1)*dimB(1)-(subB(3)-1)*lB(1))/(lB(2))+1;          
    indB = subB(1) + (subB(2)-1)*dimB(1) + (subB(3)-1)*lB(1) + (subB(4)-1)*lB(2);    
    subA(1) = subB(order(1));
    subA(2) = subB(order(2));
    subA(3) = subB(order(3));
    subA(4) = subB(order(4));
    indA = subA(1) + (subA(2)-1)*dimA(1) + (subA(3)-1)*lA(1) + (subA(4)-1)*lA(2);
    A(indA) = B(indB);
end

function A = permute5d(A, B, dimA, dimB, subA, subB, lA, lB, order, nB)

for n = 0:nB-1
    subB(1) = rem(n,dimB(1))+1;
    k = rem(n,lB(1));        
    subB(2) = (k+1-subB(1))/dimB(1)+1;  
    k = rem(n,lB(2));        
    subB(3) = (k+1-subB(1)-(subB(2)-1)*dimB(1))/(lB(1))+1;   
    k = rem(n,lB(3));            
    subB(4) = (k+1-subB(1)-(subB(2)-1)*dimB(1)-(subB(3)-1)*lB(1))/(lB(2))+1;          
    subB(5) = (n+1-subB(1)-(subB(2)-1)*dimB(1)-(subB(3)-1)*lB(1)-(subB(4)-1)*lB(2))/(lB(3))+1;      
    indB = subB(1) + (subB(2)-1)*dimB(1) + (subB(3)-1)*lB(1) + (subB(4)-1)*lB(2) + (subB(5)-1)*lB(3); 
    subA(1) = subB(order(1));
    subA(2) = subB(order(2));
    subA(3) = subB(order(3));
    subA(4) = subB(order(4));
    subA(5) = subB(order(5));
    indA = subA(1) + (subA(2)-1)*dimA(1) + (subA(3)-1)*lA(1) + (subA(4)-1)*lA(2) + (subA(5)-1)*lA(3); 
    A(indA) = B(indB);
end

function A = permute6d(A, B, dimA, dimB, subA, subB, lA, lB, order, nB)

for n = 0:nB-1
    subB(1) = rem(n,dimB(1))+1;
    k = rem(n,lB(1));        
    subB(2) = (k+1-subB(1))/dimB(1)+1;  
    k = rem(n,lB(2));        
    subB(3) = (k+1-subB(1)-(subB(2)-1)*dimB(1))/(lB(1))+1;   
    k = rem(n,lB(3));            
    subB(4) = (k+1-subB(1)-(subB(2)-1)*dimB(1)-(subB(3)-1)*lB(1))/(lB(2))+1;      
    k = rem(n,lB(4));    
    subB(5) = (k+1-subB(1)-(subB(2)-1)*dimB(1)-(subB(3)-1)*lB(1)-(subB(4)-1)*lB(2))/(lB(3))+1;      
    subB(6) = (n+1-subB(1)-(subB(2)-1)*dimB(1)-(subB(3)-1)*lB(1)-(subB(4)-1)*lB(2)-(subB(5)-1)*lB(3))/(lB(4))+1;      
    indB = subB(1) + (subB(2)-1)*dimB(1) + (subB(3)-1)*lB(1) + (subB(4)-1)*lB(2) + (subB(5)-1)*lB(3) + (subB(6)-1)*lB(4); 
    subA(1) = subB(order(1));
    subA(2) = subB(order(2));
    subA(3) = subB(order(3));
    subA(4) = subB(order(4));
    subA(5) = subB(order(5));
    subA(6) = subB(order(6));
    indA = subA(1) + (subA(2)-1)*dimA(1) + (subA(3)-1)*lA(1) + (subA(4)-1)*lA(2) + (subA(5)-1)*lA(3) + (subA(6)-1)*lA(4); 
    A(indA) = B(indB);
end

function A = permute7d(A, B, dimA, dimB, subA, subB, lA, lB, order, nB)

for n = 0:nB-1
    subB(1) = rem(n,dimB(1))+1;
    k = rem(n,lB(1));        
    subB(2) = (k+1-subB(1))/dimB(1)+1;  
    k = rem(n,lB(2));        
    subB(3) = (k+1-subB(1)-(subB(2)-1)*dimB(1))/(lB(1))+1;   
    k = rem(n,lB(3));            
    subB(4) = (k+1-subB(1)-(subB(2)-1)*dimB(1)-(subB(3)-1)*lB(1))/(lB(2))+1;      
    k = rem(n,lB(4));    
    subB(5) = (k+1-subB(1)-(subB(2)-1)*dimB(1)-(subB(3)-1)*lB(1)-(subB(4)-1)*lB(2))/(lB(3))+1;      
    k = rem(n,lB(5));    
    subB(6) = (k+1-subB(1)-(subB(2)-1)*dimB(1)-(subB(3)-1)*lB(1)-(subB(4)-1)*lB(2)-(subB(5)-1)*lB(3))/(lB(4))+1;      
    subB(7) = (n+1-subB(1)-(subB(2)-1)*dimB(1)-(subB(3)-1)*lB(1)-(subB(4)-1)*lB(2)-(subB(5)-1)*lB(3)-(subB(6)-1)*lB(4))/(lB(5))+1;      
    indB = subB(1) + (subB(2)-1)*dimB(1) + (subB(3)-1)*lB(1) + (subB(4)-1)*lB(2) + (subB(5)-1)*lB(3) + (subB(6)-1)*lB(4) + (subB(7)-1)*lB(5); 
    subA(1) = subB(order(1));
    subA(2) = subB(order(2));
    subA(3) = subB(order(3));
    subA(4) = subB(order(4));
    subA(5) = subB(order(5));
    subA(6) = subB(order(6));
    subA(7) = subB(order(7));
    indA = subA(1) + (subA(2)-1)*dimA(1) + (subA(3)-1)*lA(1) + (subA(4)-1)*lA(2) + (subA(5)-1)*lA(3) + (subA(6)-1)*lA(4) + (subA(7)-1)*lA(5); 
    A(indA) = B(indB);
end

return;

d = 4;
order = randperm(d);
dimB = randi(10,[1 d]);
B = rand(dimB);
C = permute(B,order);
[A,dimA] = cuda_permute(B, dimB, order, d);
if max(abs(A(:)-C(:)))>1e-15
    error('Array permute is wrong');
else
    disp('succsess');
end




