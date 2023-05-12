function norm = mNorm(v,i)

sizev = size(v);
n = sizev(i);

switch n
    case 2
        norm = mNorm2x2(v,i);
    case 3
        norm = mNorm3x3(v,i);
    otherwise
        norm = mNormNxN(v,i);
end

function norm = mNorm2x2(v,i)

norm = sqrt(v(i,1).^2 + v(i,2).^2);

function norm = mNorm3x3(v,i)

norm = sqrt(v(i,1).^2 + v(i,2).^2 + v(i,3).^2);

function norm = mNormNxN(v,i)

norm = sqrt(sum(v.^2,i));
