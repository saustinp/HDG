function [p,t,t2t,level] = refineall(ie, p, t, t2t, level)

if (level(ie) ~= 1),
    error('refineall error');
end
 
i1 = t(ie,1);
i2 = t(ie,2);
i3 = t(ie,3);
i4 = t(ie,4);
ie1 = t2t(ie,1);
np = size(p,1);
i5 = np+1;
i6 = np+2;
p = [p; 0.5*p(i1,1)+0.5*p(i2,1), p(i1,2)];
p = [p; 0.5*p(i1,1)+0.5*p(i2,1), p(i3,2)];

t(ie,2) = i5;
t(ie,3) = i6;
t = [t; i5, i2, i3, i6];
level = [level; level(ie)'];
np = i6;
while (ie1 > 0),
    ie = ie1;
    ie1 = t2t(ie,1);
    i1 = t(ie,1);
    i2 = t(ie,2);
    i3 = t(ie,3);
    i4 = t(ie,4);
    i6 = i5;
    i5 = np+1;
    np = i5;
    p = [p; 0.5*p(i1,1)+0.5*p(i2,1), p(i1,2)];   
    t(ie,2) = i5;
    t(ie,3) = i6;
    t = [t; i5, i2, i3, i6];
    level = [level; level(ie)'];
end

t2t = mkt2t(t,1);
