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

      
void localsolveuq(double* ae, double* fe, double* dudg, double* dudg_duh, double* M, double* C, double* E, 
        double* BD, double* F, double* GK, double* H, double* Rq, double* Ru, double* Rh, Int* ndims)
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
        
    double *Ct  = (double *) malloc( npv*npv*nd * sizeof(double) ); 
    double *Et  = (double *) malloc( ndf*npv*nd * sizeof(double) );     
    double *Bt  = (double *) malloc( npv*ncu*nd*npv*ncu * sizeof(double) ); 
    double *Dt  = (double *) malloc( npv*ncu*npv*ncu * sizeof(double) ); 
    double *Gt  = (double *) malloc( npv*ncu*nd*ncu*ndf * sizeof(double) ); 
    double *Kt  = (double *) malloc( npv*ncu*ncu*ndf * sizeof(double) ); 
    double *Ft  = (double *) malloc( ncu*ndf*npv*ncu * sizeof(double) ); 
    double *Ht  = (double *) malloc( ncu*ndf*ncu*ndf * sizeof(double) ); 
    
    Int sz[10];
    sz[0] = npv*npv;
    sz[1] = npv*ndf;
    
    // Taking transpose
    for (i=0; i<nd; i++)    
        for (j=0; j<npv; j++)                   
            for (k=0; k<npv; k++)
                Ct[k+npv*j+sz[0]*i] = C[j+npv*k+sz[0]*i];
    
    for (i=0; i<nd; i++)    
        for (j=0; j<npv; j++)                   
            for (k=0; k<ndf; k++)
                Et[k+ndf*j+sz[1]*i] = E[j+npv*k+sz[1]*i];
    
    
    sz[0] = npv*ncu;
    sz[1] = npv*ncu*npv;
    for (i=0; i<ncu; i++)    
        for (j=0; j<npv; j++)
            for (m=0; m<ncu; m++)    
                for (k=0; k<npv; k++)
                    Dt[k+npv*m+sz[0]*j+sz[1]*i] = BD[j+npv*i+sz[0]*k+sz[1]*m];        
    
    double *B;
    B = &BD[npv*ncu*npv*ncu];
        
    sz[0] = npv*ncu;
    sz[1] = npv*ncu*npv;    
    sz[2] = npv*ncu*npv*ncu;    
    for (n=0; n<nd; n++)
        for (i=0; i<ncu; i++)    
            for (j=0; j<npv; j++)            
                for (m=0; m<ncu; m++)    
                    for (k=0; k<npv; k++)
                        Bt[k+npv*m+sz[0]*j+sz[1]*i+sz[2]*n] = B[j+npv*i+sz[0]*k+sz[1]*m+sz[2]*n];        
    
    
    sz[0] = npv*ncu;
    sz[1] = npv*ncu*ncu;
    sz[2] = ncu*ndf;
    sz[3] = ncu*ndf*npv;
    for (i=0; i<ndf; i++)    
        for (j=0; j<ncu; j++)
            for (m=0; m<ncu; m++)    
                for (k=0; k<npv; k++)
                    Kt[k+npv*m+sz[0]*j+sz[1]*i] = GK[j+ncu*i+sz[2]*k+sz[3]*m];
    
    double *G;
    G = &GK[npv*ncu*ncu*ndf];
    
    sz[0] = npv*ncu;    
    sz[1] = npv*ncu*ncu;
    sz[2] = npv*ncu*ncu*ndf;
    sz[3] = ncu*ndf;
    sz[4] = ncu*ndf*npv;
    sz[5] = ncu*ndf*npv*ncu;
    for (n=0; n<nd; n++)
        for (i=0; i<ndf; i++)    
            for (j=0; j<ncu; j++)            
                for (m=0; m<ncu; m++)    
                    for (k=0; k<npv; k++)
                        Gt[k+npv*m+sz[0]*j+sz[1]*i+sz[2]*n] = G[j+ncu*i+sz[3]*k+sz[4]*m+sz[5]*n];
        
    
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
        
    sz[0] = ncu*nd;
    sz[1] = ncu*npv*ncu*nd;
    sz[2] = ncu*ncu*ndf*nd;    
    
    dpotrf(charu,&npv,M,&npv,&info);                        
    
    dpotrs(charu,&npv,&sz[0],M,&npv,Rq,&npv,&info);  
       
    dpotrs(charu,&npv,&sz[1],M,&npv,Bt,&npv,&info);  
        
    dpotrs(charu,&npv,&sz[2],M,&npv,Gt,&npv,&info);      
      
    na = ncu*ncu*ndf;
    nb = ncu*npv*ncu;    
    for (i=0; i<nd; i++) {        
        iA = i*npv*npv;
        iB = i*npv*ncu;        
        dgemm(chn, chn, &npv, &ncu, &npv, &one, &Ct[iA], &npv, &Rq[iB], &npv,
            &one, Ru, &npv);                                   
        
        iB = i*npv*ncu*npv*ncu;        
        dgemm(chn, chn, &npv, &nb, &npv, &one, &Ct[iA], &npv, &Bt[iB], &npv,
            &one, Dt, &npv);                                   
                
        iB = i*npv*ncu*ncu*ndf;
        dgemm(chn, chn, &npv, &na, &npv, &one, &Ct[iA], &npv, &Gt[iB], &npv,
            &one, Kt, &npv);                               
    }    
    
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
    
    double *dq  = (double *) malloc( npv*ncu*nd * sizeof(double) );  
    double *dlq  = (double *) malloc( npv*ncu*ncu*ndf*nd * sizeof(double) ); 
    na = ncu*ndf;
    nb = npv*ncu;
    for (i=0; i<nd; i++) {        
        iA = i*npv*ncu*npv*ncu;
        iB = i*npv*ncu;        
        dgemv(chn, &nb, &nb, &minusone, &Bt[iA], &nb, du, &inc,
            &zero, &dq[iB], &inc);                                   
        
        iB = i*npv*ncu*ncu*ndf;
        dgemm(chn, chn, &nb, &na, &nb, &minusone, &Bt[iA], &nb, dlu, &nb,
            &zero, &dlq[iB], &nb);                                                   
    }
    
    na = npv*ncu*nd;
    daxpy(&na, &one, Rq, &inc, dq, &inc);
    na = npv*ncu*ncu*ndf*nd;
    daxpy(&na, &minusone, Gt, &inc, dlq, &inc);
         
    for (i=0; i<npv*ncu; i++)
        dudg[i] = du[i];
    
    for (i=0; i<npv*ncu*nd; i++)
        dudg[npv*ncu+i] = dq[i];
            
    sz[0] = npv*ncu;
    sz[1] = npv*ncu*ncu;
    sz[2] = npv*ncu*(nd+1);
    sz[3] = npv*ncu*(nd+1)*ncu;
    for (j=0; j<ndf; j++)
        for (k=0; k<ncu; k++)
            for (m=0; m<ncu; m++)
                for (n=0; n<npv; n++)
                    dudg_duh[n+npv*m+sz[2]*k+sz[3]*j] = dlu[n+npv*m+sz[0]*k+sz[1]*j];
    
    
    sz[0] = npv*ncu;
    sz[1] = npv*ncu*(nd+1);
    sz[2] = npv*ncu*(nd+1)*ncu;
    sz[3] = npv*ncu*ncu;
    sz[4] = npv*ncu*ncu*ndf;
    sz[5] = npv*ndf;    
    for (j=0; j<ndf; j++)
        for (k=0; k<ncu; k++)
            for (i=0; i<nd; i++)
                for (m=0; m<ncu; m++)
                    for (n=0; n<npv; n++) 
                        dudg_duh[n+npv*m+sz[0]*(i+1)+sz[1]*k+sz[2]*j] 
                          = dlq[n+npv*m+sz[0]*k+sz[3]*j+sz[4]*i]; 
    
    
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
    
    double *fet  = (double *) malloc( ndf*ncu * sizeof(double) );  
    double *aet  = (double *) malloc( ndf*ncu*ncu*ndf * sizeof(double) ); 

    Int colsB = ncu*ncu*ndf;
    dgemm(chn, chn, &ndf, &ncu, &npv, &one, &Et[0], &ndf, &dq[0], 
                &npv, &zero, &fet[0], &ndf);    
    dgemm(chn, chn, &ndf, &colsB, &npv, &one, &Et[0], &ndf, &dlq[0], 
                &npv, &zero, &aet[0], &ndf);        
    for (i=1; i<nd; i++) {     
        iA = i*ndf*npv;        
        iB = i*npv*ncu;        
        dgemm(chn, chn, &ndf, &ncu, &npv, &one, &Et[iA], &ndf, &dq[iB], 
                    &npv, &one, &fet[0], &ndf);    
        
        iB = i*npv*ncu*ncu*ndf;        
        dgemm(chn, chn, &ndf, &colsB, &npv, &one, &Et[iA], &ndf, &dlq[iB], 
                    &npv, &one, &aet[0], &ndf);    
    }        
    
    for (j=0; j<ndf; j++)
        for (i=0; i<ncu; i++)
            fe[i+ncu*j] = fe[i+ncu*j] - fet[j+ndf*i];

    sz[0] = ncu*ndf;
    sz[1] = ncu*ndf*ncu;
    for (n=0; n<ndf; n++)
        for (m=0; m<ncu; m++)    
            for (j=0; j<ndf; j++)
                for (i=0; i<ncu; i++)
                    ae[i+ncu*j+sz[0]*m+sz[1]*n] = ae[i+ncu*j+sz[0]*m+sz[1]*n] + aet[j+ndf*i+sz[0]*m+sz[1]*n];    
    
    free(Et); free(Ct); free(Dt); free(Bt); 
    free(Gt); free(Kt); free(Ft); free(Ht);        
    free(Rt); free(ipiv); free(dq); free(dlq);     
    free(fet); free(aet);
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
    
    
    Int iA, iB, iC, iD, iE, iF, iG, iH, iI, iJ, iK, iL;        
    for (i=0; i<ne; i++) {        
        iA = i*npv*nc;   
        iB = i*npv*nc*ncu*ndf;        
        iC = i*npv*npv;        
        iD = i*npv*npv*nd;
        iE = i*npv*ndf*nd;                                        
        iF = i*npv*ncu*npv*nc;                        
        iG = i*npv*ncu*ncu*ndf;       
        iH = i*npv*ncu*nd;        
        iI = i*npv*ncu;                
        iJ = i*ncu*ndf*ncu*ndf;
        iK = i*ncu*ndf;
        iL = i*ncu*ndf*npv*nc;
        
        localsolveuq(&AE[iJ], &FE[iK], &DUDG[iA], &DUDG_DUH[iB], &M[iC], &C[iD], &E[iE], 
                &BD[iF], &F[iG], &GK[iL], &H[iJ], &Rq[iH], &Ru[iI], &Rh[iK], ndims);                   
        
    }            
}
    
