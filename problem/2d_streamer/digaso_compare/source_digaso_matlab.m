function [s,s_udg] = source_digaso_matlab(s, s_udg, pg, udg, param, time, ng, nc, ncu, nd, ncd)

    % Note that the arrays are pre-flattened before passing into the
    % function

    % Include the header for the C code
    coder.cinclude('source_streamer2d.h');

    % Evaluate the C function - the source function returns void
    coder.ceval('source_streamer2d', coder.ref(s), coder.ref(s_udg), coder.ref(pg), coder.ref(udg), coder.ref(param), int32(time), int32(ng), int32(nc), int32(ncu), int32(nd), int32(ncd)); 
    
    end
