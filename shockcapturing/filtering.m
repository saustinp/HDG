function Uout = filtering(Uin, nodelist, weight)

[np,ncu] = size(Uin);
Uout = zeros(np,ncu);

for i = 1:np    
    % list of nodes around node i       
    node = nodelist{i}; 
    % field at those nodes
    un = Uin(node,:);
    % weight at those nodes
    w = weight{i};
    % weighted sum of the field at those nodes
    Uout(i,:) = sum(un.*w,1)/sum(w);        
end


