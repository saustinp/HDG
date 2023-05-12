function [p,t,t2t,level] = refine(ie, p, t, t2t, level)

i1 = t(ie,1);
i2 = t(ie,2);
i3 = t(ie,3);
i4 = t(ie,4);
ie1 = t2t(ie,1);
np = size(p,1);
i5 = np+1;
i6 = np+2;
i7 = np+3;
i8 = np+4;
p = [p; 2*p(i1,1)/3+p(i2,1)/3, p(i1,2)];
p = [p; p(i1,1)/3+2*p(i2,1)/3, p(i1,2)];
p = [p; p(i1,1)/3+2*p(i2,1)/3, p(i1,2)/2+p(i3,2)/2];
p = [p; 2*p(i1,1)/3+p(i2,1)/3, p(i1,2)/2+p(i3,2)/2];

t(ie,2) = i5;
t(ie,3) = i8;
t = [t; i8, i7, i3, i4];
t = [t; i5, i6, i7, i8];
t = [t; i6, i2, i3, i7];
level = [level; level(ie)*[1,1,1]'];
np = i8;
while (ie1 > 0),
    ie = ie1;
    ie1 = t2t(ie,1);
    i1 = t(ie,1);
    i2 = t(ie,2);
    i3 = t(ie,3);
    i4 = t(ie,4);
    i7 = i6;
    i8 = i5;
    i5 = np+1;
    i6 = np+2;
    np = i6;
    p = [p; 2*p(i1,1)/3+p(i2,1)/3, p(i1,2)];
    p = [p; p(i1,1)/3+2*p(i2,1)/3, p(i1,2)];
    
    t(ie,2) = i5;
    t(ie,3) = i8;
    t = [t; i5, i6, i7, i8];
    t = [t; i6, i2, i3, i7];
    level = [level; level(ie)*[1,1]'];
end

t2t = mkt2t(t,1);
