#ifndef __FHATDRIVER
#define __FHATDRIVER

//#include "fluxDriver.cpp"
#include "FM/stabilizationTensor.cpp"

// Written by: C. Nguyen & P. Fernandez

void UDG_2_UDG4fh(double *UDG_aux, double * UDG, double * UH, Int fhatExpression, int numPoints, int ncu, int ncq, int nc)
{
    Int inc = 1, len;

    if (fhatExpression == 0) {
        /* FHdotN = dot(F(UH,Q),n) + tau*(U-UH). These are the usual trace fluxes in HDG */
        
        len = numPoints * ncu;
        DCOPY(&len, &UH[0], &inc, &UDG_aux[0], &inc);

        len = numPoints * ncq;
        DCOPY(&len, &UDG[ncu * numPoints], &inc, &UDG_aux[ncu * numPoints], &inc);
        /* UDG_aux: numPoints / nc */
        /* UDG: numPoints / nc */
        /* UH: numPoints / nch */
    }
    else if (fhatExpression == 1) {
        /* FHdotN = dot(F(U,Q),n) + tau*(U-UH). These fluxes are used in interfaces with LAMBDA playing the role of UH. */

        len = numPoints * nc;
        DCOPY(&len, &UDG[0], &inc, &UDG_aux[0], &inc);
        /* UDG_aux: numPoints / nc */
        /* UDG: numPoints / nc */
    }
    else if (fhatExpression == 2) {
        /* FHdotN = dot(F(U,Q),n) */
        
        len = numPoints * nc;
        DCOPY(&len, &UDG[0], &inc, &UDG_aux[0], &inc);
        /* UDG_aux: numPoints / nc */
        /* UDG: numPoints / nc */
    }
    else {
        printf("fhatExpression = %d has invalid value in function UDG_2_UDG4fh\n", fhatExpression);
        printf("Execution will be terminated\n");
        exit(-1);
    }
}



void fh2fhn(double* fhn, double* fhn_UDG, double* fhn_UH, double* fh, double* fh_UDG, double* UDG, double* UH,
            double* tau, double* tau_UH, double* nl, int numPoints, int ncu, int ncq, int nc, int nd, Int fhatExpression, int computeJacobian)
{
    // TODO: In part 2 of this function, introduce an alternative implementation for constant tau and symmetric stabilization tensors.

    /* Note: Trying to use BLAS here to perform the dot product with the normal vector would be slower than a for loop. */

    int i, j, k, l, n;
    int nm, nk, nn, ifhn;
    int sz2, sz3;

    // Part 1: Contribution due to fh*nl
    sz2 = numPoints * ncu;
    for (j=0; j<ncu; j++)
        for (i=0; i<numPoints; i++) {
            ifhn = j * numPoints + i;
            fhn[ifhn] = fh[0*sz2 + j * numPoints + i] * nl[i + 0 * numPoints];
            for (l = 1; l < nd; l++)
                fhn[ifhn] += fh[l*sz2 + j * numPoints + i] * nl[i + l * numPoints];
        }
    // fhn: numPoints / ncu
    // fh: numPoints / ncu / nd
    // nl: numPoints / nd

    if (computeJacobian == 1) {
        int ifhn_UDG, ifh_UDG, ifhn_UH;

        if (fhatExpression == 0) {
            // fhn = dot(F(UH,Q),n) + tau*(U-UH). These are the usual trace fluxes in HDG

            for (k = 0; k < ncu; k++)
                for (j = 0; j < ncu; j++)
                    for (i = 0; i < numPoints; i++)
                        fhn_UDG[i + j * numPoints + k * sz2] = 0.0;

            sz3 = sz2 * nd;
            for (k = 0; k < ncq; k++)
                for (j = 0; j < ncu; j++)
                    for (i = 0; i < numPoints; i++) {
                        ifhn_UDG = i + j * numPoints + (ncu + k) * sz2;
                        ifh_UDG = i + j * numPoints + (ncu + k) * sz3;
                        fhn_UDG[ifhn_UDG] = fh_UDG[ifh_UDG + 0 * sz2] * nl[i + 0 * numPoints];
                        for (l = 1; l < nd; l++)
                            fhn_UDG[ifhn_UDG] += fh_UDG[ifh_UDG + l * sz2] * nl[i + l * numPoints];
                    }
            // fhn_UDG: numPoints / ncu / nc
            // fh_UDG: numPoints / ncu / nd / nc
            // nl: numPoints / nd

            for (k = 0; k < ncu; k++)
                for (j = 0; j < ncu; j++)
                    for (i = 0; i < numPoints; i++) {
                        ifhn_UH = i + j * numPoints + k * sz2;
                        ifh_UDG = i + j * numPoints + k * sz3;
                        fhn_UH[ifhn_UH] = fh_UDG[ifh_UDG + 0 * sz2] * nl[i + 0 * numPoints];
                        for (l = 1; l < nd; l++) {
                            fhn_UH[ifhn_UH] += fh_UDG[ifh_UDG + l * sz2] * nl[i + l * numPoints];
                        }
                    }
            // fhn_UH: numPoints / ncu / nch
            // fh_UDG: numPoints / ncu / nd / nch
            // nl: numPoints / nd
        }
        else if (fhatExpression == 1 || fhatExpression == 2) {
            // fhatExpression == 1: FHdotN = dot(F(U,Q),n) + tau*(U-UH). These fluxes are used in interfaces with LAMBDA playing the role of UH
            // fhatExpression == 2: FHdotN = dot(F(U,Q),n)
            
            sz3 = sz2 * nd;
            for (k = 0; k < nc; k++)
                for (j = 0; j < ncu; j++)
                    for (i = 0; i < numPoints; i++) {
                        ifhn_UDG = i + j * numPoints + k * sz2;
                        ifh_UDG = i + j * numPoints + k * sz3;
                        fhn_UDG[ifhn_UDG] = fh_UDG[ifh_UDG + 0 * sz2] * nl[i + 0 * numPoints];
                        for (l = 1; l < nd; l++) {
                            fhn_UDG[ifhn_UDG] += fh_UDG[ifh_UDG + l * sz2] * nl[i + l * numPoints];
                        }
                    }
            // fhn_UDG: numPoints / ncu / nc
            // fh_UDG: numPoints / ncu / nd / nc
            // nl: numPoints / nd

            int sz_fhn_UH = numPoints * ncu * ncu;
            for (i = 0; i < sz_fhn_UH; i++)
                fhn_UH[i] = 0.0;
            // fhn_UH: numPoints / ncu / nch
            // fh_UH: numPoints / ncu / nd / nch
            // nl: numPoints / nd
        }
        else {
            printf("fhatExpression = %d has invalid value in function fh2fhn\n", fhatExpression);
            printf("Execution will be terminated\n");
            exit(-1);
        }
    }

    // Part 2: Contribution due to tau*(U-UH)
    if (fhatExpression == 0 || fhatExpression == 1) {
        sz2 = ncu * numPoints;
        sz3 = ncu * ncu * numPoints;
        for (k=0; k<ncu; k++)
            for (j=0; j<ncu; j++)
                for (i=0; i<numPoints; i++) {
                    nm = k * sz2 + j * numPoints + i;
                    nk = k * numPoints + i;
                    fhn[j * numPoints + i] += tau[nm] * (UDG[nk] - UH[nk]);

                    if (computeJacobian == 1) {
                        fhn_UDG[nm] += tau[nm];
                        fhn_UH[nm] -= tau[nm];
                        for (n = 0; n < ncu; n++) {
                            nk = k * sz3 + n * sz2 + j * numPoints + i;
                            nn = n * numPoints + i;
                            fhn_UH[nm] += tau_UH[nk] * (UDG[nn] - UH[nn]);
                        }
                    }
                }
        // fhn: numPoints / ncu
        // fhn_UDG: numPoints / ncu / nc
        // fhn_UH: numPoints / ncu / nch
        // tau: numPoints / ncu / ncu
        // tau_UH: numPoints / ncu / ncu / nch
        // UDG: numPoints / nc
        // UH: numPoints / nch
    }
}



void fhatDriver(double * fhn, double * fhn_UDG, double * fhn_UH, double * pg, double * UDG, double * UH, double * UDGref,
                double * UHref, double * NL, double *avField_p1CG, meshstruct &mesh, masterstruct &master,
                appstruct &app, double * param, double time, Int* ndims, Int fhatExpression, int computeJacobian,
                int numPoints)
{
    int nfe, ncu, nch, ncq, nc, nd, ncd;
    nd  = (int) ndims[0];
    ncd = (int) ndims[1];
    nfe = (int) ndims[2];
    nc  = (int) ndims[19];
    ncu = (int) ndims[20];
    ncq = (int) ndims[21];
    nch = (int) ndims[23];
    
    double *UDG_aux = new double[numPoints * nc];
    double *UDG_aux_ref = new double[numPoints * nc];
    double *fh = new double[numPoints * ncu * nd];
    double *fh_UDG = new double[numPoints * ncu * nd * nc];
    double *tau = new double[numPoints * ncu * ncu];
    double *tau_UH = new double[numPoints * ncu * ncu * ncu];
    
    UDG_2_UDG4fh(UDG_aux, UDG, UH, fhatExpression, numPoints, ncu, ncq, nc);
    UDG_2_UDG4fh(UDG_aux_ref, UDGref, UHref, fhatExpression, numPoints, ncu, ncq, nc);
    
    fluxDriver(fh, fh_UDG, pg, UDG_aux, UDG_aux_ref, avField_p1CG, mesh, master, app, param, time, ndims, computeJacobian, numPoints);
    
//     getStabilizationTensor(tau, tau_UH, UH, pg, NL, app, param, numPoints, nch, nd, computeJacobian);
    getStabilizationTensor(tau, tau_UH, UHref, pg, NL, app, param, numPoints, nch, nd, computeJacobian);
    
    fh2fhn(fhn, fhn_UDG, fhn_UH, fh, fh_UDG, UDG, UH, tau, tau_UH, NL, numPoints, ncu, ncq, nc, nd, fhatExpression, computeJacobian);

    delete[] UDG_aux; delete[] UDG_aux_ref; delete[] fh; delete[] fh_UDG; delete[] tau; delete[] tau_UH;
}

#endif
