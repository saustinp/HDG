
#include "../fbou.cpp"

// Written by: C. Nguyen & P. Fernandez

void fbou_nsNEW(double *fh, double *fh_u, double *fh_uh, double *pg, double *udg, double *udg_ref, double *uhg,
                double *uhg_ref, double *nl, double *avg_p1CG, double *ui, meshstruct &mesh,
                masterstruct &master, appstruct &app, double *param, double time, Int ib, Int* ndims,
                int computeJacobian, int numPoints)
{
    fbou(fh, fh_u, fh_uh, pg, udg, udg_ref, uhg, uhg_ref, nl, avg_p1CG, ui, mesh, master, app, param, time, ib, ndims, computeJacobian, numPoints);
}
