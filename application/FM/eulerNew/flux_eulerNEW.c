
#include "flux_euler2dNEW.c"
#include "flux_euler3dNEW.c"

// Written by: C. Nguyen & P. Fernandez

void flux_eulerNEW(double *f, double *f_udg, double *pg, double *udg, double *udg_ref, double *avField_p1CG, meshstruct &mesh, masterstruct &master, appstruct &app, double *param,
                double time, int ng, int nc, int ncu, int nd, int ncd, Int* ndims, int computeJacobian)
{
    switch (nd) {
        case 2:
            flux_euler2dNEW(f, f_udg, pg, udg, app, param, time, ng, nc, ncu, nd, ncd, computeJacobian);
            // Add artificial viscosity if necessary
            break;
        case 3:
            flux_euler3dNEW(f, f_udg, pg, udg, app, param, time, ng, nc, ncu, nd, ncd, computeJacobian);
            // Add artificial viscosity if necessary
            break;
        default:
            exit(-1);
            break;
    }
}
