function uhat = hdg_linearsystem(Ae, Re, elcon, neumann, uc, c)

[npfe, nt] = size(Re);

il = zeros(npfe,npfe,nt);
jl = zeros(npfe,npfe,nt);
for i=1:nt
    con = repmat((elcon(:,i)'-1),1,1)+ones(1,npfe);
    con = reshape(con,npfe,1);
    il(:,:,i) = repmat(con ,1,npfe);
    jl(:,:,i) = repmat(con',npfe,1);        
end
H = sparse(reshape(il,(npfe)^2*nt,1),reshape(jl,(npfe)^2*nt,1),reshape(Ae,(npfe)^2*nt,1));
R = sparse(reshape(il(:,1,:),(npfe)*nt,1),ones((npfe)*nt,1),reshape(Re,(npfe)*nt,1));                    

if neumann==1
  R(end+1) = -c;
  H(end+1,:) = uc;  
  H(:,end+1) = [uc(:); 4*max(abs(uc(:)))];
  ut = H\R;
  uhat = ut(1:end-1);    
else  
  uhat = H\R;
end
