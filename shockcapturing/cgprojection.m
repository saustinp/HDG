function [u, udg] = cgprojection(ucg, A, B, cgelcon, cgent2dgent, rowent2elem) 

[npe,ne] = size(cgelcon);
nc = size(ucg,2);
nqe = size(B,1);

% convert CG field to DG field 
udg = reshape(ucg(cgelcon, :), [npe ne*nc]);
tm = udg;

% convert nodal DG field to modal DG field 
udg = reshape(A\udg, [npe ne*nc]);

if (nqe>npe) % prolongation operator
    udg(npe+1:nqe,:) = 0;    
else % restriction operator
    udg(nqe+1:npe,:) = [];    
end

%[tm(:,51) (A*udg(1:npe,51))]

% convert modal DG field to nodal DG field
udg = permute(reshape(B*udg, [nqe ne nc]),[1 3 2]);

% convert DG field to CG field
u = udg2ucg(udg, cgent2dgent, rowent2elem);

end


