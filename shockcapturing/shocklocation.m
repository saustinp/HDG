function [px, ux] = shocklocation(mesh, s, idx, x, y, nref)

[m, n] = size(x);
sx = dgfieldatx(mesh,s,idx,[x(:) y(:)],nref);
sx = reshape(sx, [m, n]);

nd = mesh.nd;
px = zeros(m,nd);
ux = zeros(m,1);
for i = 1:m
    [smax, imax] = max(sx(i,:));
    ux(i) = smax;
    px(i,:) = [x(i,imax) y(i, imax)];
end

