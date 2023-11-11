function [fh,fh_udg, fh_uh] = fhat_digaso_matlab(fh, fh_udg, fh_uh, pg, udg, uh, nl, param, time, ng, nc, nch, nd, ncd)

    % Note that the arrays are pre-flattened before passing into the
    % function

    % Include the header for the C code
    coder.cinclude('fhat_streamer2d.h');

    % Evaluate the C function - the source function returns void
    coder.ceval('fhat_streamer2d', coder.ref(fh), coder.ref(fh_udg), coder.ref(fh_uh), coder.ref(pg), coder.ref(udg), coder.ref(uh), coder.ref(nl), coder.ref(param), int32(time), int32(ng), int32(nc), int32(nch), int32(nd), int32(ncd));
   
    end
