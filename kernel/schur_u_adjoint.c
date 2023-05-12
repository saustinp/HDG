#include "mex.h"
#include "blas.h"
#include "lapack.h"

typedef ptrdiff_t     Int;

void print1darray(double* a, Int m)
{
    Int i;
    for (i=0; i<m; i++) 
      mexPrintf("%g  ",a[i]);
    mexPrintf("\n");
    mexPrintf("\n");
}

void print2darray(double* a, Int m, Int n)
{
    Int i, j;
    for (i=0; i<m; i++) {
        for (j=0; j<n; j++) 
            mexPrintf("%g  ",a[j*m+i]);
        mexPrintf("\n");
    }    
    mexPrintf("\n");
}

void print3darray(double* a, Int m, Int n, Int p)
{
    Int i, j, k;
    
    for (k=0; k<p; k++) {
        for (i=0; i<m; i++) {
            for (j=0; j<n; j++) 
                mexPrintf("%g  ",a[k*n*m+j*m+i]);
            mexPrintf("\n");
        }
        mexPrintf("\n");
    }
    mexPrintf("\n");
}

      
void localsolve(double* ae, double* fe, double* dudg, double* dudg_duh,
        double* BD, double* F, double* GK, double* H, double* Ru, double* Rh, Int* ndims)
{         
    Int inc = 1, i, j, m, n, k, is, ks, iA, iB, iC, na, nb, info;
    char *chn = "N", *charu = "U";
    double one = 1.0, zero = 0.0, minusone = -1.0, fac, tm;        
    
    /* Get dimensions */
    Int npv, ncu, nc, ndf;    
    nc  = ndims[0]; 
    ncu = ndims[1];      
    npv = ndims[2]; 
    ndf = ndims[3];        
        
        
    double *Dt  = (double *) malloc( npv*ncu*npv*ncu * sizeof(double) );     
    double *Kt  = (double *) malloc( npv*ncu*ncu*ndf * sizeof(double) ); 
    double *Ft  = (double *) malloc( ncu*ndf*npv*ncu * sizeof(double) ); 
    double *Ht  = (double *) malloc( ncu*ndf*ncu*ndf * sizeof(double) ); 
    
    Int sz[10];
    
    // Taking transpose    
    sz[0] = npv*ncu;
    sz[1] = npv*ncu*npv;
    for (i=0; i<ncu; i++)    
        for (j=0; j<npv; j++)
            for (m=0; m<ncu; m++)    
                for (k=0; k<npv; k++)
                    Dt[k+npv*m+sz[0]*j+sz[1]*i] = BD[j+npv*i+sz[0]*k+sz[1]*m];        
    
    
    sz[0] = npv*ncu;
    sz[1] = npv*ncu*ncu;
    sz[2] = ncu*ndf;
    sz[3] = ncu*ndf*npv;
    for (i=0; i<ndf; i++)    
        for (j=0; j<ncu; j++)
            for (m=0; m<ncu; m++)    
                for (k=0; k<npv; k++)
                    Kt[k+npv*m+sz[0]*j+sz[1]*i] = GK[j+ncu*i+sz[2]*k+sz[3]*m];
        
    sz[0] = ncu*ndf;
    sz[1] = ncu*ndf*npv;
    sz[2] = npv*ncu;
    sz[3] = npv*ncu*ncu;
    for (i=0; i<ncu; i++)    
        for (j=0; j<npv; j++)
            for (m=0; m<ndf; m++)    
                for (k=0; k<ncu; k++)
                    Ft[k+ncu*m+sz[0]*j+sz[1]*i] = F[j+npv*i+sz[2]*k+sz[3]*m];
    
    
    sz[0] = ncu*ndf;
    sz[1] = ncu*ndf*ncu;    
    for (i=0; i<ndf; i++)    
        for (j=0; j<ncu; j++)
            for (m=0; m<ndf; m++)    
                for (k=0; k<ncu; k++)
                    Ht[k+ncu*m+sz[0]*j+sz[1]*i] = H[j+ncu*i+sz[0]*k+sz[1]*m];
        
          
    double *Rt  = (double *) malloc( npv*ncu*(1+ndf*ncu) * sizeof(double) ); 
    for (i=0; i<npv*ncu; i++)
        Rt[i] = Ru[i];
    
    for (j=0; j<ndf*ncu; j++)
        for (i=0; i<npv*ncu; i++)        
            Rt[(j+1)*npv*ncu+i] = -Kt[j*npv*ncu+i];
    
    Int *ipiv =  (Int *) malloc( npv*ncu * sizeof(Int) );    
    na = npv*ncu;
    nb = 1+ndf*ncu;
    dgesv(&na, &nb, Dt, &na, ipiv, Rt, &na, &info);
             
    double *du, *dlu;    
    du  = &Rt[0];
    dlu = &Rt[npv*ncu];
    
    na = npv*ncu;
    dcopy(&na,du,&inc,dudg,&inc);    
                        
    na = npv*ncu*ncu*ndf;
    dcopy(&na,dlu,&inc,dudg_duh,&inc);
        
    na = ncu*ndf;
    dcopy(&na,Rh,&inc,fe,&inc);
    
    na = ncu*ndf*ncu*ndf;
    dcopy(&na,Ht,&inc,ae,&inc);
                
    na = ncu*ndf;
    nb = npv*ncu;
    dgemv(chn, &na, &nb, &minusone, Ft, &na, du, 
                    &inc, &one, fe, &inc);               
    dgemm(chn, chn, &na, &na, &nb, &one, Ft, &na, dlu, 
                    &nb, &one, ae, &na);    
        
    free(Dt); free(Kt); free(Ft); free(Ht);        
    free(Rt); free(ipiv);
}

void mexFunction(int nlhs,mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
            
    const mwSize *szA, *szB, *szC;    
    Int npv, ncu, nc, ndf, ne, ndim;    
    Int i, inc = 1;
    char *chn = "N";
    double one = 1.0, minusone = -1.0;  
        
    double* BD = mxGetPr(prhs[0]);
    double* F  = mxGetPr(prhs[1]);
    double* GK = mxGetPr(prhs[2]);
    double* H  = mxGetPr(prhs[3]);
    double* Ru = mxGetPr(prhs[4]);       
    double* Rh = mxGetPr(prhs[5]);       
        
    szA = mxGetDimensions(prhs[0]);                
    npv = szA[0];
    ncu = szA[1];
    nc  = szA[3];
    ne  = szA[4];     
    ndim = mxGetNumberOfDimensions(prhs[0]);
    if (ndim == 4) {
        ne = 1;
    }
        
    szB = mxGetDimensions(prhs[3]);                
    ndf = szB[1];
        
    Int ndims[5];        
    ndims[0] = nc; 
    ndims[1] = ncu;      
    ndims[2] = npv; 
    ndims[3] = ndf;        
    
    mwSize szS[3]; szS[0]=npv; szS[1]=nc; szS[2]=ne; 
    plhs[0]=mxCreateNumericArray(3,szS,mxDOUBLE_CLASS,mxREAL);         
    
    mwSize szT[5]; szT[0]=npv; szT[1]=nc; szT[2]=ncu; szT[3]=ndf; szT[4]=ne; 
    plhs[1]=mxCreateNumericArray(5,szT,mxDOUBLE_CLASS,mxREAL);         
    
    szT[0]=ncu; szT[1]=ndf; szT[2]=ncu; szT[3]=ndf; szT[4]=ne; 
    plhs[2]=mxCreateNumericArray(5,szT,mxDOUBLE_CLASS,mxREAL);         
    
    szS[0]=ncu; szS[1]=ndf; szS[2]=ne; 
    plhs[3]=mxCreateNumericArray(3,szS,mxDOUBLE_CLASS,mxREAL);         
    
    double* DUDG = mxGetPr(plhs[0]);
    double* DUDG_DUH = mxGetPr(plhs[1]);    
    double* AE = mxGetPr(plhs[2]);
    double* FE = mxGetPr(plhs[3]);    
    
    Int iA, iB, iC, iD, iE, iF, iG, iH, iI, iJ, iK, iL;        
    for (i=0; i<ne; i++) {        
        iA = i*npv*nc;   
        iB = i*npv*nc*ncu*ndf;        
        iF = i*npv*ncu*npv*nc;                        
        iG = i*npv*ncu*ncu*ndf;       
        iI = i*npv*ncu;                
        iJ = i*ncu*ndf*ncu*ndf;
        iK = i*ncu*ndf;
        iL = i*ncu*ndf*npv*nc;        
        localsolve(&AE[iJ], &FE[iK], &DUDG[iA], &DUDG_DUH[iB], 
                   &BD[iF], &F[iG], &GK[iL], &H[iJ], &Ru[iI], &Rh[iK], ndims);                   
        
    }            
}
    
