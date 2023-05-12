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

      
void localsolve(double* dudg, double* dudg_duh, double* BD, double* F, double* Ru, Int* ndims)
{         
    Int inc = 1, i, j, na, nb, info;
    char *chn = "N", *charu = "U";    
    
    /* Get dimensions */
    Int npv, ncu, nc, ndf;       
    nc  = ndims[0]; 
    ncu = ndims[1];      
    npv = ndims[2]; 
    ndf = ndims[3];        
                        
    double *Rt  = (double *) malloc( npv*ncu*(1+ndf*ncu) * sizeof(double) ); 
    
    for (i=0; i<npv*ncu; i++)
        Rt[i] = Ru[i];
    
    for (j=0; j<ndf*ncu; j++)
        for (i=0; i<npv*ncu; i++)        
            Rt[(j+1)*npv*ncu+i] = -F[j*npv*ncu+i];
    
    Int *ipiv =  (Int *) malloc( npv*ncu * sizeof(Int) );
    
    na = npv*ncu;
    nb = 1+ndf*ncu;
    dgesv(&na, &nb, BD, &na, ipiv, Rt, &na, &info);
             
    double *du, *dlu;
    du  = &Rt[0];
    dlu = &Rt[npv*ncu];           
        
    na = npv*ncu;
    dcopy(&na,du,&inc,dudg,&inc);    
                        
    na = npv*ncu*ncu*ndf;
    dcopy(&na,dlu,&inc,dudg_duh,&inc);
                             
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
        
    Int N;
    N = ncu*ndf*ncu*ndf*ne;
    dcopy(&N,&H[0],&inc,&AE[0],&inc);
    N = ncu*ndf*ne;
    dcopy(&N,&Rh[0],&inc,&FE[0],&inc);
    
    Int iA, iB, iC, iD, iE, iF, iG, iH, iI, rowsA, colsA;    
    rowsA = ncu*ndf;
    colsA = npv*nc;
    for (i=0; i<ne; i++) {        
        iA = i*npv*nc;   
        iB = i*npv*nc*ndf*ncu;        
        iF = i*npv*npv*ncu*nc;                        
        iG = i*npv*ncu*ndf*ncu;               
        iI = i*npv*ncu;                
         
        localsolve(&DUDG[iA], &DUDG_DUH[iB], &BD[iF], &F[iG], &Ru[iI], ndims);                   
                
        iA = i*rowsA*colsA;
        iB = i*colsA*rowsA;
        iC = i*rowsA*rowsA;
        dgemm(chn, chn, &rowsA, &rowsA, &colsA, &one, &GK[iA], &rowsA, &DUDG_DUH[iB], 
                    &colsA, &one, &AE[iC], &rowsA);
        
        iB = i*colsA;
        iC = i*rowsA;
        dgemv(chn, &rowsA, &colsA, &minusone, &GK[iA], &rowsA, &DUDG[iB], 
                    &inc, &one, &FE[iC], &inc);                           
    }            
}
    
