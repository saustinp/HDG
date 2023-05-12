function elemR = elementlistcg(mesh,cgnodes,R)
% for each cg node
%   find all elements within a distance R from the node

[ne,nv] = size(mesh.t);
[np,nd] = size(cgnodes);

% element centers
p = reshape(mesh.p(mesh.t',:),[nv ne nd]);
xm = reshape(mean(p,1),[ne nd]);

R2 = R*R;
elemR = cell(np,1);
if nd==2
    for i = 1:np        
        s2 = (cgnodes(i,1)-xm(:,1)).^2 + (cgnodes(i,2)-xm(:,2)).^2;    
        elemR{i} = find(s2 <= R2);
    end
elseif nd==3
    for i = 1:ne
        s2 = (cgnodes(i,1)-xm(:,1)).^2 + (cgnodes(i,2)-xm(:,2)).^2 + (cgnodes(i,3)-xm(:,3)).^2;    
        elemR{i} = find(s2 <= R2);
    end    
end

