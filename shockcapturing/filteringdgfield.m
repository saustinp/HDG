function ucg = filteringdgfield(udg, dgnodes, elemlist, R, wR)

[npe,ncu,ne] = size(udg);
ucg = zeros(npe,ncu,ne);
for i = 1:ne
    xm = dgnodes(:,:,i);
    elem = elemlist{i};    
    nn = length(elem);
    for j = 1:npe   
        tm = 0;        
        wm = 0;        
        for k = 1:nn
            un = udg(:,:,elem(k));            
            xn = dgnodes(:,:,elem(k));            
            d = sqrt((xm(j,1)-xn(:,1)).^2 + (xm(j,2)-xn(:,2)).^2);                            
            idx = (d <= R);
            if sum(idx) > 0
                w  = 1 - d(idx)*((1 - wR)/R);       
                wm = wm + sum(w);
                tm = tm + sum(un(idx).*w,1);
            end
        end        
        ucg(j,:,i) = tm/wm;
    end
end

end
   
% function w = weight(s,wR,R)
%     % s = 0 -> w = 1 
%     % s = R -> w = wR
%     w  = 1 - s*((1 - wR)/R);
% end

