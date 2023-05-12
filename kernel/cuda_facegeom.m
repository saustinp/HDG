function [pg, nlg, jac] = cuda_facegeom(pg, nlg, jac, dpg, pn, shapt, nf, ng, np, nd, ncx)
% VOLGEOM computes physical nodes, Jacobian determinant and inverse at Gauss points  

%   [pg, Xx, jac] = volgeom(shapgeomvt,dgnodes)
%
%    SHAPGEOMVT :  Shape functions and derivatives at Gauss points
%    DGNODES    :  Geometry DG nodes 
%    PG         :  Physical nodes at Gauss points 
%    Xx         :  Minus the inverse of the Jacobian mapping times the determinant
%    jac        :  Determinant of the Jacobian mapping 

pg = cuda_node2gauss(pg, pn, shapt, ng, np, nf*ncx);
for i = 1:nd-1
    dpg(:,:,i) = cuda_node2gauss(dpg(:,:,i), pn, shapt(:,:,i+1), ng, np, nf*nd);
end
dpg = reshape(dpg,[ng*nf nd nd-1]);

switch nd
    case 1
        for i = 1:ng*nf
            jac(i) = 1;
            nlg(i) = 1;
        end
    case 2
        for i = 1:ng*nf
            jac(i) = sqrt(dpg(i,1).^2+dpg(i,2).^2);
            nlg(i,1)= dpg(i,2)/jac(i);            
            nlg(i,2) = -dpg(i,1)/jac(i);
        end
    case 3
        for i = 1:ng*nf
            nlg(i,1) = dpg(i,2,1).*dpg(i,3,2) - dpg(i,3,1).*dpg(i,2,2);
            nlg(i,2) = dpg(i,3,1).*dpg(i,1,2) - dpg(i,1,1).*dpg(i,3,2);
            nlg(i,3) = dpg(i,1,1).*dpg(i,2,2) - dpg(i,2,1).*dpg(i,1,2);
            jac(i) = sqrt(nlg(i,1).^2+nlg(i,2).^2+nlg(i,3).^2);
            nlg(i,1) = nlg(i,1)/jac(i);
            nlg(i,2) = nlg(i,2)/jac(i);
            nlg(i,3) = nlg(i,3)/jac(i);            
        end
    otherwise
        error('Dimension is not implemented');
end
dpg = reshape(dpg,[ng*nf nd*(nd-1)]);
