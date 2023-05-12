function [fhg, fhg_udgp, fhg_udgm] = cuda_fhat(ib, udg1, udg2, uhg, uinf, xdg1, nlg1, jac1, param, time, ng, ncu, nd, nc)
%TODO: the provided cuda_fhat uses ldgfhat and ldgfbou, which do not have Jacobian terms. 
%      The functions euler_fhat and euler_fbou do have jacobian terms, but I think that
%      euler_fbou is specific to HDG and needs to be written in a ldg form. 
    if ib==0
    %     fhg = ldgfhat(nlg1, xdg1, udg1, udg2, uhg, param, time);
        [fhg, fhg_udgp, fhg_udgm] = euler_fhat(nlg1,xdg1,udg1,udg2,param,time);
    elseif ib>0
        fhg = ldgfbou(ib, uinf(1:ncu), nlg1, xdg1, udg1, uhg, param, time);
    end
    
    for i = 1:ncu
        fhg(:,i) = fhg(:,i).*jac1(:);
    end
    
    for i = 1:ncu
        for j = 1:ncu
            fhg_udgp(:,i,j) = fhg_udgp(:,i,j).*jac1(:);
            fhg_udgm(:,i,j) = fhg_udgm(:,i,j).*jac1(:);
        end
    end
    
    