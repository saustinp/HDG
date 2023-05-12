function elemR = elementlist(mesh,R)

nd = mesh.nd;
[ne,nv] = size(mesh.t);

% element centers
p = reshape(mesh.p(mesh.t',:),[nv ne nd]);
xm = reshape(mean(p,1),[ne nd]);

R2 = R*R;
elemR = cell(ne,1);
if nd==2
    for i = 1:ne
        s2 = (xm(i,1)-xm(:,1)).^2 + (xm(i,2)-xm(:,2)).^2;    
        elemR{i} = find(s2 <= R2);
    end
elseif nd==3
    for i = 1:ne
        s2 = (xm(i,1)-xm(:,1)).^2 + (xm(i,2)-xm(:,2)).^2 + (xm(i,3)-xm(:,3)).^2;    
        elemR{i} = find(s2 <= R2);
    end    
end

