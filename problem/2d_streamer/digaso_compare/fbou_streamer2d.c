#include <math.h>
#include "fhat_streamer2d.c"

void fbou_streamer2d(double *fh, double *fh_u, double *fh_uh, 
          double *pg, double *udg, double *uhg, double *nl,
          double *ui, double *param, double time, int ib,
          int ng, int nc, int nch, int nd, int ncd)
{            
    
    int    i;                   
    double kappa   = param[0];
    double tau = param[1];   
    double r;       // Mod for axisymmetric
    
    if (ib==1) { // Electrode: species set to homogeneous nuemann and potential set to dirichlet
        
        // Zero all arrays
        for (i=0; i<ng*nch; i++)
            fh[i] = 0.0;
        for (i=0; i<ng*nch*nc; i++)
            fh_u[i] = 0.0;
        for (i=0; i<ng*nch*nch; i++)
            fh_uh[i] = 0.0;

        fhat_streamer2d(fh, fh_u, fh_uh, pg, udg, uhg, nl, param, time, ng, nch, nc, nd, ncd);

        // fh(:,3) = 0
        for (i=0; i<ng; i++)
            fh[2*ng + i] = 0.0;

        // fh_udg(:,3,:) = 0
        for (int k=0; k<nc; k++){
            for (int j=0; j<nch; j++){
                for (int i=0; i<ng; i++) {
                    if (j == 2)
                        fh_u[ng*nch*k + ng*j + i] = 0.0;
                }
            }
        }

        // fh_uh(:,3,:) = 0
        for (int k=0; k<nch; k++){
            for (int j=0; j<nch; j++){
                for (int i=0; i<ng; i++) {
                    if (j == 2)
                        fh_uh[ng*nch*k + ng*j + i] = 0.0;
                }
            }
        }

        for (i=0; i<ng; i++) {
            
            // Potential dirichlet
            // l_ref = param[0];
            // E_ref = param[2];
            // phi0 = param[4];
            // phi0_tilde = phi0/(E_ref*l_ref);

            fh[2*ng+i] = r*tau*(62.5-uhg[2*ng+i]);     // NOTE: phi0_tilde hardcoded for now
            fh_uh[2*ng*nch+2*ng+i] = -r*tau;
        }
  
    }                                               
    else if (ib==2) { // Right farfield -- species + potential all have homogeneous neumann
        fhat_streamer2d(fh, fh_u, fh_uh, pg, udg, uhg, nl, param, time, ng, nch, nc, nd, ncd);    

    }                                     
    else if (ib==3) { // Symmetry boundary: extrapolate m=u or u=uhat
        // Zero all arrays
        for (i=0; i<ng*nch; i++)
            fh[i] = 0.0;
        for (i=0; i<ng*nch*nc; i++)
            fh_u[i] = 0.0;
        for (i=0; i<ng*nch*nch; i++)
            fh_uh[i] = 0.0;

        for (i=0; i<ng; i++) {

            // fh(:,1), fh(:,2), fh(:,3)
            fh[0*ng+i] = r*tau*(udg[0*ng+i]-uhg[0*ng+i]);        // Check that uhg is correct here?
            fh[1*ng+i] = r*tau*(udg[1*ng+i]-uhg[1*ng+i]);
            fh[2*ng+i] = r*tau*(udg[2*ng+i]-uhg[2*ng+i]);
            
            // fh_udg(:,1,1), fh_udg(:,2,2), fh_udg(:,3,3)
            fh_u[0*ng*nch+0*ng+i] = r*tau;
            fh_u[1*ng*nch+1*ng+i] = r*tau;
            fh_u[2*ng*nch+2*ng+i] = r*tau;

            // fh_uh(:,1,1), fh_uh(:,2,2), fh_uh(:,3,3)
            fh_uh[0*ng*nch+0*ng+i] = -r*tau;
            fh_uh[1*ng*nch+1*ng+i] = -r*tau;
            fh_uh[2*ng*nch+2*ng+i] = -r*tau;
        }       
    }                           
    else {                        
        // printf("This error is in %s on line %d\n",__FILE__, __LINE__);
        // printf("Boundary condition %d is not implemented yet.", ib);            
        // exit(-1);                                    
    }                
}