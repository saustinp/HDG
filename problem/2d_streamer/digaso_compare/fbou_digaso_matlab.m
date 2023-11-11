function [fh,fh_u, fh_uh] = fbou_digaso_matlab(fh, fh_u, fh_uh, pg, udg, uhg, nl, ui, param, time, ib, ng, nc, nch, nd, ncd)
    % Note that the arrays are pre-flattened before passing into the
    % function

    % Include the header for the C code
    coder.cinclude('fbou_streamer2d.h');

    % Evaluate the C function - the source function returns void
    coder.ceval('fbou_streamer2d', coder.ref(fh), coder.ref(fh_u), coder.ref(fh_uh), coder.ref(pg), coder.ref(udg), coder.ref(uhg), coder.ref(nl), coder.ref(ui), coder.ref(param), int32(time), int32(ib), int32(ng), int32(nc), int32(nch), int32(nd), int32(ncd));
   
    end
