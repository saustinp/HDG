function [pl, pr] = boundingbox(px)

[~, nd] = size(px);

if nd==1
    pl = min(px);
    pr = max(px);    
elseif nd==2
    x1 = min(px(:,1));
    x2 = max(px(:,1));
    y1 = min(px(:,2));
    y2 = max(px(:,2));
    pl = [x1 y1];
    pr = [x2 y2];
elseif nd==3
    x1 = min(px(:,1));
    x2 = max(px(:,1));
    y1 = min(px(:,2));
    y2 = max(px(:,2));
    z1 = min(px(:,3));
    z2 = max(px(:,3));
    pl = [x1 y1 z1];
    pr = [x2 y2 z2];    
end

