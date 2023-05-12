function [pp, pl, pr, pa] = shockcells(px, d, n, m)

[np,nd] = size(px);
[pl, pr] = boundingbox(px);
pa = zeros(m, nd, n+1);

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
            if i==1
                r = linspace(rmin, rmax, m);
            else
                r = linspace(min(rmin,pp(i-1,3)), max(rmax, pp(i-1,4)), m);                
            end
            pa(:,:,i) = [s(i)*ones(m,1) r]; 
        end
        pa(:,:,n+1) = [s(end)*ones(m,1) r]; 
    else
        c = 1;
        for i = 1:n
            ind = (px(:,d) >= s(i)) & (px(:,d) <= s(i+1));
            rmin = min(px(ind,c));
            rmax = max(px(ind,c));
            pp(i,:) = [rmin rmax s(i) s(i+1)];                         
            if i==1                
                r = linspace(rmin, rmax, m);
            else
                r = linspace(min(rmin,pp(i-1,1)), max(rmax, pp(i-1,2)), m);                
            end
            pa(:,:,i) = [r(:) s(i)*ones(m,1)]; 
        end                
        pa(:,:,n+1) = [r(:) s(end)*ones(m,1)]; 
    end
end

