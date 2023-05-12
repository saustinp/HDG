function Uout = filteringcgfield(Uin, cgnodes, nodelist, R, wR)

[np,ncu] = size(Uin);
Uout = zeros(np,ncu);

for i = 1:np
    xm = cgnodes(i,:);  % node i
    node = nodelist{i}; % list of nodes around node i       
    xn = cgnodes(node,:); 
    % distance from node i to those nodes around it
    d = sqrt((xm(1)-xn(:,1)).^2 + (xm(2)-xn(:,2)).^2);    
    % distance-based weights
    w  = 1 - d*((1 - wR)/R);     
    % field at those nodes
    un = Uin(node,:);
    % weighted sum of the field at those nodes
    Uout(i,:) = sum(un.*w,1)/sum(w);        
end

