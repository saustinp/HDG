function [out] = apply_weight_to_vector(W, X, mesh, master)
% Used for applying a block matrix (Mass or Mass inverse) to a vector: WX
% Used for weighted LSPG or for POD. 
% NOTE/TODO: if a computational bottleneck, can be implemented in C code 

if isequal(length(W), 1) %scalar case
    out = W*X;
else

    npe = master.npe;
%     ncu = pde.ncu;
    ne = mesh.ne;
    ncu = length(X) / (ne*npe); %ncu can change depending on whether we are applying POD to system or components
    u_shape = [npe ncu ne];
    u_flat = [npe*ncu*ne 1];
    out = zeros(u_shape);
    
    if isequal(size(X),u_flat)
        X = reshape(X, u_shape);
        reshape_flag = 1;
    end
    
    for ie = 1:ne
        out(:,:,ie) = W(:,:,ie) * X(:,:,ie);
    end
    
    if reshape_flag
        out = reshape(out, u_flat);
    end
end
end