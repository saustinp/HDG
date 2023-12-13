function [f,f_udg] = flux_digaso_matlab(f, f_udg, pg, udg, param, time, ng, nc, ncu, nd, ncd)

    % Note that the arrays are pre-flattened before passing into the
    % function

    % Include the header for the C code
    coder.cinclude('flux.h');

    % Evaluate the C function - the source function returns void
    coder.ceval('flux', coder.ref(f), coder.ref(f_udg), coder.ref(pg), coder.ref(udg), coder.ref(param), int32(time), int32(ng), int32(nc), int32(ncu), int32(nd), int32(ncd));
   
    end
