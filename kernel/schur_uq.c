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

      
void localsolveuq(double* dudg, double* dudg_duh, double* M, double* C, double* E, 
        double* BD, double* F, double* Rq, double* Ru, Int* ndims)
{         
    Int inc = 1, i, j, m, n, k, is, ks, iA, iB, iC, na, nb, info;
    char *chn = "N", *charu = "U";
    double one = 1.0, zero = 0.0, minusone = -1.0, fac, tm;        
    
    /* Get dimensions */
    Int nd, npv, ncu, nc, ndf, nd1;
    nd  = ndims[0];    
    nc  = ndims[1]; 
    ncu = ndims[2];      
    npv = ndims[3]; 
    ndf = ndims[4];    
    nd1 = nd+1;
                
    double *B, *D;
    
    Int sz[10];
    sz[0] = ncu*nd;
    sz[1] = npv*nd;
    sz[2] = ndf*nd;    
    
    D = &BD[0];
    B = &BD[npv*ncu*npv*ncu];
        
    /* compute the inverse of M using Cholesky decomposition */
    dpotrf(charu,&npv,M,&npv,&info);                        
    
    /* compute inv(M)*Rq and store it in Rq */        
    dpotrs(charu,&npv,&sz[0],M,&npv,Rq,&npv,&info);  
    
    /* compute inv(M)*Ce and store it in C */        
    dpotrs(charu,&npv,&sz[1],M,&npv,C,&npv,&info);  
    
    /* compute inv(M)*Me and store it in E */    
    dpotrs(charu,&npv,&sz[2],M,&npv,E,&npv,&info);               
    
    double *BMiR  = (double *) malloc( npv*ncu * sizeof(double) ); 
    double *BMiC  = (double *) malloc( npv*ncu*npv*ncu * sizeof(double) ); 
    double *BMiE  = (double *) malloc( npv*ncu*ndf*ncu * sizeof(double) ); 
    
    na = npv*ncu;
    nb = npv*ncu*nd;        
    dgemv(chn, &na, &nb, &one, B, &na, Rq, 
                    &inc, &zero, BMiR, &inc);                        
    
    double *zd  = (double *) malloc( nd * sizeof(double) ); 
    zd[0]=0.0;
    for (i=1; i<nd; i++)
        zd[i] = 1.0;
    
    for (j=0; j<ncu; j++)
        for (i=0; i<nd; i++) {        
            iA = i*npv*ncu*npv*ncu+j*npv*ncu*npv;
            iB = i*npv*npv;
            iC = j*npv*ncu*npv;
            na = npv*ncu;        
            dgemm(chn, chn, &na, &npv, &npv, &one, &B[iA], &na, &C[iB], &npv,
                &zd[i], &BMiC[iC], &na);                                   

            iB = i*npv*ndf;   
            iC = j*npv*ncu*ndf;
            dgemm(chn, chn, &na, &ndf, &npv, &one, &B[iA], &na, &E[iB], &npv,
                &zd[i], &BMiE[iC], &na);                                   
        }    
                        
    na = npv*ncu*npv*ncu;
    daxpy(&na, &one, BMiC, &inc, D, &inc);
    
    na = npv*ncu;
    daxpy(&na, &minusone, BMiR, &inc, Ru, &inc);
        
    sz[0] = ncu*npv;
    sz[1] = ncu*ncu*npv;
    sz[3] = npv*ncu*ndf;
    for (j=0; j<ndf; j++)
        for (k=0; k<ncu; k++)
            for (m=0; m<ncu; m++)
                for (n=0; n<npv; n++)                     
                    F[n+npv*m+sz[0]*k+sz[1]*j] = BMiE[n+npv*m+sz[0]*j+sz[3]*k] - F[n+npv*m+sz[0]*k+sz[1]*j];                    
    
    double *Rt  = (double *) malloc( npv*ncu*(1+ndf*ncu) * sizeof(double) ); 
    for (i=0; i<npv*ncu; i++)
        Rt[i] = Ru[i];
    
    for (j=0; j<ndf*ncu; j++)
        for (i=0; i<npv*ncu; i++)        
            Rt[(j+1)*npv*ncu+i] = F[j*npv*ncu+i];
           
    Int *ipiv =  (Int *) malloc( npv*ncu * sizeof(Int) );
    
    na = npv*ncu;
    nb = 1+ndf*ncu;
    dgesv(&na, &nb, D, &na, ipiv, Rt, &na, &info);
             
    double *du, *dlu;
    du  = &Rt[0];
    dlu = &Rt[npv*ncu];           
    
    
    double *MiCdu  = (double *) malloc( npv*ncu*nd * sizeof(double) ); 
    double *MiCdlu = (double *) malloc( npv*ncu*ndf*ncu*nd * sizeof(double) ); 
    
    for (i=0; i<nd; i++) {        
        iA = i*npv*npv;
        iB = i*npv*ncu;
        dgemm(chn, chn, &npv, &ncu, &npv, &one, &C[iA], &npv, du, &npv,
            &zero, &MiCdu[iB], &npv);                                   
        
        na = ncu*ndf*ncu;
        iB = i*npv*ncu*ncu*ndf;
        dgemm(chn, chn, &npv, &na, &npv, &one, &C[iA], &npv, dlu, &npv,
            &zero, &MiCdlu[iB], &npv);                                                   
    }
    
    for (i=0; i<npv*ncu; i++)
        dudg[i] = du[i];
    
    for (i=0; i<npv*ncu*nd; i++)
        dudg[npv*ncu+i] = Rq[i] + MiCdu[i];
                    
    sz[2] = npv*ncu*nd1;
    sz[3] = npv*ncu*nd1*ncu;
    for (j=0; j<ndf; j++)
        for (k=0; k<ncu; k++)
            for (m=0; m<ncu; m++)
                for (n=0; n<npv; n++)
                    dudg_duh[n+npv*m+sz[2]*k+sz[3]*j] = dlu[n+npv*m+sz[0]*k+sz[1]*j];
    
    sz[0] = npv*ncu;
    sz[1] = npv*ncu*nd1;
    sz[2] = npv*ncu*nd1*ncu;
    sz[3] = npv*ncu*ncu;
    sz[4] = npv*ncu*ncu*ndf;
    sz[5] = npv*ndf;    
    for (j=0; j<ndf; j++)
        for (k=0; k<ncu; k++)
            for (i=0; i<nd; i++)
                for (m=0; m<ncu; m++)
                    for (n=0; n<npv; n++) {
                        fac = (m == k) ? 1.0 : 0.0;
                        dudg_duh[n+npv*m+sz[0]*(i+1)+sz[1]*k+sz[2]*j] 
                          = MiCdlu[n+npv*m+sz[0]*k+sz[3]*j+sz[4]*i] - fac*E[n+npv*j+sz[5]*i]; 
                    }
                     
    free(zd); free(BMiR); free(BMiC); free(BMiE);
    free(Rt); free(ipiv);
    free(MiCdu); free(MiCdlu);
}

void mexFunction(int nlhs,mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
            
    const mwSize *szA, *szB, *szC;    
    Int nd, npv, ncu, nc, ndf, ne, ndim;    
    Int i, inc = 1;
    char *chn = "N";
    double one = 1.0, minusone = -1.0;  
        
    double* M  = mxGetPr(prhs[0]);                       
    double* C  = mxGetPr(prhs[1]);       
    double* E  = mxGetPr(prhs[2]);
    double* BD = mxGetPr(prhs[3]);
    double* F  = mxGetPr(prhs[4]);
    double* GK = mxGetPr(prhs[5]);
    double* H  = mxGetPr(prhs[6]);
    double* Rq = mxGetPr(prhs[7]);
    double* Ru = mxGetPr(prhs[8]);       
    double* Rh = mxGetPr(prhs[9]);       
        
    szA = mxGetDimensions(prhs[0]);                
    npv = szA[0];    
    ne  = szA[2]; 
    
    ndim = mxGetNumberOfDimensions(prhs[0]);
    if (ndim == 2) {
        ne = 1;
    }
        
    szB = mxGetDimensions(prhs[2]);                
    ndf = szB[1];
    nd  = szB[2];    
    
    szC = mxGetDimensions(prhs[3]);                
    ncu = szC[1];
    nc  = szC[3];
        
    Int ndims[5];    
    ndims[0] = nd;    
    ndims[1] = nc; 
    ndims[2] = ncu;      
    ndims[3] = npv; 
    ndims[4] = ndf;        
    
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
        iC = i*npv*npv;        
        iD = i*npv*npv*nd;
        iE = i*npv*ndf*nd;                                        
        iF = i*npv*npv*ncu*nc;                        
        iG = i*npv*ncu*ndf*ncu;       
        iH = i*npv*ncu*nd;        
        iI = i*npv*ncu;                
        localsolveuq(&DUDG[iA], &DUDG_DUH[iB], &M[iC], &C[iD], &E[iE], 
                &BD[iF], &F[iG], &Rq[iH], &Ru[iI], ndims);                   
                
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
    
