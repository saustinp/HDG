function fg = cuda_flux(xdgg, udgg, param, time, ng, ncu, nd, nc)

    fg = euler_flux( xdgg, udgg, param, time);
    