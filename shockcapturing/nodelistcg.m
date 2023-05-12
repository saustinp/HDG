function [nodeR, weightR] = nodelistcg(mesh,cgnodes,cgelcon,R,wR)
% for each cg node
%   find all CG nodes within a distance R from the node

[ne,nv] = size(mesh.t);
[np,nd] = size(cgnodes);

% element centers
p = reshape(mesh.p(mesh.t',:),[nv ne nd]);
pm = reshape(mean(p,1),[ne nd]);

R2 = R*R;
nodeR = cell(np,1);
weightR = cell(np,1);
if nd==2
    for i = 1:np        
        % get node i
        x = cgnodes(i,1);
        y = cgnodes(i,2);
        % distance from node i to element centers
        s2 = (x-pm(:,1)).^2 + (y-pm(:,2)).^2;    
        % find all elements within a distance R from (x, y)
        elem = s2 <= R2; 
        % get all CG nodes associated with those elements
        node = cgelcon(:,elem);        
        node = unique(node(:));   
        n = length(node);
        idx = zeros(n,1);
        % find all elements within a distance R from (x, y)
        for j = 1:n % loop over all CG nodes found
            xj = cgnodes(node(j),1);
            yj = cgnodes(node(j),2);
            s2 = (x-xj)^2 + (y-yj)^2;    
            if (s2 <= R2)
                idx(j) = 1;
            end
        end        
        node = node(idx==1);
        nodeR{i} = node;        
        xn = cgnodes(node,:); 
        % distance from node i to those nodes around it
        d = sqrt((x-xn(:,1)).^2 + (y-xn(:,2)).^2);    
        % distance-based weights
        weightR{i}   = 1 - d*((1 - wR)/R);             
    end
elseif nd==3
    for i = 1:ne
        % get node i
        x = cgnodes(i,1);
        y = cgnodes(i,2);
        z = cgnodes(i,3);        
        s2 = (x-pm(:,1)).^2 + (y-pm(:,2)).^2 + (z-pm(:,3)).^2;    
        elem = s2 <= R2;
        % get all CG nodes associated with those elements
        node = cgelcon(:,elem);        
        node = unique(node(:));   
        n = length(node);
        idx = zeros(n,1);
        % find all elements within a distance R from (x, y)
        for j = 1:n % loop over all CG nodes found
            xj = cgnodes(node(j),1);
            yj = cgnodes(node(j),2);
            zj = cgnodes(node(j),3);
            s2 = (x-xj)^2 + (y-yj)^2 + (z-zj)^2;    
            if (s2 <= R2)
                idx(j) = 1;
            end
        end
        node = node(idx==1);
        nodeR{i} = node;        
        xn = cgnodes(node,:); 
        % distance from node i to those nodes around it
        d = sqrt((x-xn(:,1)).^2 + (y-xn(:,2)).^2 + (z-xn(:,3)).^2);    
        % distance-based weights
        weightR{i}   = 1 - d*((1 - wR)/R);             
    end    
end

