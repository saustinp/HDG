function [pg, Xx, jac] = cuda_elemgeom(pg, Xx, jac, pn, Jg, shapt, ne, ng, np, nd, ncx)
% VOLGEOM computes physical nodes, Jacobian determinant and inverse at Gauss points  

%   [pg, Xx, jac] = volgeom(shapgeomvt,dgnodes)
%
%    SHAPGEOMVT :  Shape functions and derivatives at Gauss points
%    DGNODES    :  Geometry DG nodes 
%    PG         :  Physical nodes at Gauss points 
%    Xx         :  Minus the inverse of the Jacobian mapping times the determinant
%    jac        :  Determinant of the Jacobian mapping 

Jg = reshape(Jg, [ng*ne nd nd]);
Xx = reshape(Xx, [ng*ne nd nd]);

pg = cuda_node2gauss(pg, pn, shapt, ng, np, ne*ncx);
for i = 1:nd
    Jg(:,:,i) = cuda_node2gauss(Jg(:,:,i), pn, shapt(:,:,i+1), ng, np, ne*nd);
end

switch nd
    case 1
        for i = 1:ng*ne
            jac(i) = Jg(i);
            Xx(i) = 1;
        end
    case 2
        for i = 1:ng*ne
            jac(i) = Jg(i,1,1).*Jg(i,2,2) - Jg(i,2,1).*Jg(i,1,2);  
            Xx(i,1,1) = Jg(i,2,2);
            Xx(i,2,1) = -Jg(i,1,2);
            Xx(i,1,2) = -Jg(i,2,1);
            Xx(i,2,2) = Jg(i,1,1);            
        end
    case 3
        for i = 1:ng*ne
            jac(i) = Jg(i,1,1).*Jg(i,2,2).*Jg(i,3,3) - Jg(i,1,1).*Jg(i,2,3).*Jg(i,3,2)+ ...
                     Jg(i,1,2).*Jg(i,2,3).*Jg(i,3,1) - Jg(i,1,2).*Jg(i,2,1).*Jg(i,3,3)+ ...
                     Jg(i,1,3).*Jg(i,2,1).*Jg(i,3,2) - Jg(i,1,3).*Jg(i,2,2).*Jg(i,3,1);            
            Xx(i,1,1) = Jg(i,2,2).*Jg(i,3,3) - Jg(i,3,2).*Jg(i,2,3);
            Xx(i,2,1) = Jg(i,3,2).*Jg(i,1,3) - Jg(i,1,2).*Jg(i,3,3);
            Xx(i,3,1) = Jg(i,1,2).*Jg(i,2,3) - Jg(i,2,2).*Jg(i,1,3);
            Xx(i,1,2) = Jg(i,3,1).*Jg(i,2,3) - Jg(i,2,1).*Jg(i,3,3);
            Xx(i,2,2) = Jg(i,1,1).*Jg(i,3,3) - Jg(i,3,1).*Jg(i,1,3);
            Xx(i,3,2) = Jg(i,2,1).*Jg(i,1,3) - Jg(i,1,1).*Jg(i,2,3);
            Xx(i,1,3) = Jg(i,2,1).*Jg(i,3,2) - Jg(i,3,1).*Jg(i,2,2);
            Xx(i,2,3) = Jg(i,3,1).*Jg(i,1,2) - Jg(i,1,1).*Jg(i,3,2);
            Xx(i,3,3) = Jg(i,1,1).*Jg(i,2,2) - Jg(i,2,1).*Jg(i,1,2);                            
        end                
    otherwise
        error('Dimension is not implemented');
end

Jg = reshape(Jg, [ng*ne nd*nd]);
Xx = reshape(Xx, [ng*ne nd*nd]);

