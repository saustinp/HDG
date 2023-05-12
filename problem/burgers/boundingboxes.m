function [pp, pl, pr, pa] = boundingboxes(px, d, n, k)

[~, nd] = size(px);
[pl, pr] = boundingbox(px);
pa = zeros(k(1)*k(2),nd,n);

if nd==2   
    smin = pl(d);
    smax = pr(d);
    s = linspace(smin, smax, n+1);
    pp = zeros(n,4);
    if d==1
        c = 2;
        for i = 1:n
            ind = (px(:,d) >= s(i)) & (px(:,d) <= s(i+1));
            rmin = min(px(ind,c));
            rmax = max(px(ind,c));
            pp(i,:) = [s(i) s(i+1) rmin rmax];            
        end
    else
        c = 1;
        for i = 1:n
            ind = (px(:,d) >= s(i)) & (px(:,d) <= s(i+1));
            rmin = min(px(ind,c));
            rmax = max(px(ind,c));
            pp(i,:) = [rmin rmax s(i) s(i+1)];                         
        end                
    end    
    for i = 1:n
        x = linspace(pp(i,1),pp(i,2),k(1));
        y = linspace(pp(i,3),pp(i,4),k(2));
        [X,Y] = ndgrid(x,y);
        pa(:,:,i) = [X(:) Y(:)];
    end
end
pa = reshape(pa,[k nd n]);

