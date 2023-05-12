function [px, ux] = shockcurve(mesh, s, s0, nref)

A = refinementmatrix(mesh,nref);

npe = size(s,1);
nc = size(s,2);
nd = mesh.nd;
npl = size(A,1);

% find all elements containing shock
sm = squeeze(max(s,[],1));
idx = find(sm > s0);
n = length(idx);

u = s(:,:,idx);
p = mesh.dgnodes(:,1:nd,idx);

uref = reshape(A*reshape(u,npe,nc*n),[npl nc n]);
pref = reshape(A*reshape(p,npe,nd*n),[npl nd n]);

un = zeros(n,1);
pn = zeros(n,nd);
for i = 1:n
    ui = uref(:,1,i);
    pi = pref(:,:,i);
    [umax, imax] = max(ui);
    pn(i,:) = pi(imax,:);
    un(i) = umax;
end

% remove close nodes
k = 1;
px = 0*pn;
ux = 0*un;
[~,idx] = max(un);
ux(1) = un(idx);
px(1,:) = pn(idx,:);
for i = 1:n
    in = 0;
    for j = 1:k
        if  norm(pn(i,:) - px(j,:)) < 1e-6
            in = 1;
            break;
        end
    end
    if in == 0        
        k = k + 1;
        px(k,:) = pn(i,:);
        ux(k) = un(i);
    end
end
px = px(1:(k-1),:);
ux = ux(1:(k-1));

idx = ones(k-1,1);
for i = 1:(k-1)
    for j = 1:(k-1)
        if (i ~= j)
            if norm(px(i,:) - px(j,:)) < 0.04
                if ux(i) < ux(j)
                    idx(i) = 0;
                else
                    idx(j) = 0;
                end
            end
        end
    end
end
px = px(idx==1,:);
ux = ux(idx==1);






