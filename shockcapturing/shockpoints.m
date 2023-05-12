function [px, sx, idx] = shockpoints(mesh, s, s0, s1, nref)

A = refinementmatrix(mesh,nref);

npe = size(s,1);
nd = mesh.nd;
npl = size(A,1);

% find all elements containing shock
sm = squeeze(mean(s,1));
idx = find(sm > s0);
n = length(idx);

u = s(:,idx);
p = permute(mesh.dgnodes(:,1:nd,idx),[1 3 2]);

sref = reshape(A*reshape(u,npe,n),[npl*n 1]);
pref = reshape(A*reshape(p,npe,n*nd),[npl*n nd]);

ind = sref > s1;
sx = sref(ind);
px = pref(ind,:);

