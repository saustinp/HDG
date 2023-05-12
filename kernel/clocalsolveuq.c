#include "mex.h"
#include "blas.h"
#include "lapack.h"

typedef ptrdiff_t     Int;


void localsolveuq(double* dudg, double* dudg_duh, double* Ru, double* BD, double* Rlu, 
        double* MiC, double* MiE, double* udg, double* uh, double* sh, double fc_q, Int* ndims)
{         
    Int inc = 1, i, j, m, n, k, is, ks, iA, iB, iC, na, nb;
    char *chn = "N";
    double one = 1.0, zero = 0.0, minusone = -1.0, fac, tm;        
    
    /* Get dimensions */
    Int nd, npv, nch, nc, ndf, nd1;
    nd  = ndims[0];    
    nc  = ndims[1]; 
    nch = ndims[2];      
    npv = ndims[3]; 
    ndf = ndims[4];    
    nd1 = nd+1;
                
    double *u, *q, *s, *B, *D;
    
    Int sz[8];
    sz[0] = npv*npv;
    sz[1] = nch*sz[0];
    sz[2] = nch*npv;
    sz[3] = npv*sz[2];
    
    u = &udg[0];
    q = &udg[npv*nch];
    s = &sh[npv*nch];
    D = &BD[0];
    B = &BD[npv*nch*npv*nch];
    
    fac = 1/fc_q;
    na = npv*npv*nd;
    dscal(&na,&fac,MiC,&inc);
        
    na = npv*ndf*nd;
    dscal(&na,&fac,MiE,&inc);
    
    double *MiCu  = (double *) malloc( npv*nch*nd * sizeof(double) ); 
    double *MiEu  = (double *) malloc( npv*nch*nd * sizeof(double) ); 
    
    for (i=0; i<nd; i++) {        
        iA = i*npv*npv;
        iB = i*npv*nch;
        dgemm(chn, chn, &npv, &nch, &npv, &one, &MiC[iA], &npv, u, &npv,
            &zero, &MiCu[iB], &npv);                                   
        
        iA = i*npv*ndf;        
        dgemm(chn, chn, &npv, &nch, &ndf, &one, &MiE[iA], &npv, uh, &ndf,
            &zero, &MiEu[iB], &npv);                                   
    }
            
    double *MiR  = (double *) malloc( npv*nch*nd * sizeof(double) ); 
    double *BMiR  = (double *) malloc( npv*nch * sizeof(double) ); 
    double *BMiC  = (double *) malloc( npv*nch*npv*nch * sizeof(double) ); 
    double *BMiE  = (double *) malloc( npv*nch*ndf*nch * sizeof(double) ); 
    
    for (i=0; i<npv*nch*nd; i++)
        MiR[i] = s[i]*fac + MiEu[i] - MiCu[i] - q[i]; 
        
    na = npv*nch;
    nb = npv*nch*nd;        
    dgemv(chn, &na, &nb, &minusone, B, &na, MiR, 
                    &inc, &zero, BMiR, &inc);            
        
    double *zd  = (double *) malloc( nd * sizeof(double) ); 
    zd[0]=0.0;
    for (i=1; i<nd; i++)
        zd[i] = 1.0;
    
    for (j=0; j<nch; j++)
        for (i=0; i<nd; i++) {        
            iA = i*npv*nch*npv*nch+j*npv*nch*npv;
            iB = i*npv*npv;
            iC = j*npv*nch*npv;
            na = npv*nch;        
            dgemm(chn, chn, &na, &npv, &npv, &minusone, &B[iA], &na, &MiC[iB], &npv,
                &zd[i], &BMiC[iC], &na);                                   

            iB = i*npv*ndf;   
            iC = j*ndf*nch*npv;
            dgemm(chn, chn, &na, &ndf, &npv, &minusone, &B[iA], &na, &MiE[iB], &npv,
                &zd[i], &BMiE[iC], &na);                                   
        }    
                    
    na = npv*nch*npv*nch;
    daxpy(&na, &one, BMiC, &inc, D, &inc);
    
    na = npv*nch;
    daxpy(&na, &one, BMiR, &inc, Ru, &inc);
    
    na = npv*nch*ndf*nch;
    daxpy(&na, &one, BMiE, &inc, Rlu, &inc);   
    
    double *Rt  = (double *) malloc( npv*nch*(1+ndf*nch) * sizeof(double) ); 
    for (i=0; i<npv*nch; i++)
        Rt[i] = Ru[i];
    
    for (j=0; j<ndf*nch; j++)
        for (i=0; i<npv*nch; i++)        
            Rt[(j+1)*npv*nch+i] = Rlu[j*npv*nch+i];
           
    Int info;
    Int *ipiv =  (Int *) malloc( npv*nch * sizeof(Int) );
    
    na = npv*nch;
    nb = 1+ndf*nch;
    dgesv(&na, &nb, D, &na, ipiv, Rt, &na, &info);
        
    double *du, *dlu;
    du  = &Rt[0];
    dlu = &Rt[npv*nch];
    
    double *MiCdu  = (double *) malloc( npv*nch*nd * sizeof(double) ); 
    double *MiCdlu = (double *) malloc( npv*nch*ndf*nch*nd * sizeof(double) ); 
    
    for (i=0; i<nd; i++) {        
        iA = i*npv*npv;
        iB = i*npv*nch;
        dgemm(chn, chn, &npv, &nch, &npv, &one, &MiC[iA], &npv, du, &npv,
            &zero, &MiCdu[iB], &npv);                                   
        
        na = nch*ndf*nch;
        iB = i*npv*nch*ndf*nch;
        dgemm(chn, chn, &npv, &na, &npv, &one, &MiC[iA], &npv, dlu, &npv,
            &zero, &MiCdlu[iB], &npv);                                                   
    }
    
    for (i=0; i<npv*nch; i++)
        dudg[i] = du[i];
    
    for (i=0; i<npv*nch*nd; i++)
        dudg[npv*nch+i] = MiR[i] - MiCdu[i];
        
    sz[0] = nch*npv;
    sz[1] = nd1*sz[0];
    sz[2] = ndf*sz[1];
    sz[3] = ndf*sz[0];    
    sz[4] = nch*sz[3];
    sz[5] = ndf*npv;
    for (j=0; j<nch; j++)
        for (k=0; k<ndf; k++)
            for (m=0; m<nch; m++)
                for (n=0; n<npv; n++)
                    dudg_duh[j*sz[2]+k*sz[1]+m*npv+n] = dlu[j*sz[3]+k*sz[0]+m*npv+n];
    
    for (i=0; i<nd; i++)
        for (j=0; j<nch; j++)
            for (k=0; k<ndf; k++)
                for (m=0; m<nch; m++)
                    for (n=0; n<npv; n++) {
                        fac = (m == j) ? 1.0 : 0.0;
                        dudg_duh[j*sz[2]+k*sz[1]+(i+1)*sz[0]+m*npv+n] 
                                = fac*MiE[i*sz[5]+k*npv+n] - MiCdlu[i*sz[4]+j*sz[3]+k*sz[0]+m*npv+n];                
                    }
               
    free(Rt); free(zd); free(ipiv);
    free(MiCu); free(MiEu); free(MiR);
    free(BMiR); free(BMiC); free(BMiE);     
    free(MiCdu); free(MiCdlu); 
        
}

void mexFunction(int nlhs,mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
            
    const mwSize *szA, *szB, *szC;    
    Int nd, npv, nch, nc, ndf, ne, ndim;    
    Int i;
    
    double* BD = mxGetPr(prhs[0]);                       
    double* Ru = mxGetPr(prhs[1]);       
    double* Rlu = mxGetPr(prhs[2]);
    double* MiC = mxGetPr(prhs[3]);
    double* MiE = mxGetPr(prhs[4]);
    double* UDG = mxGetPr(prhs[5]);
    double* UH = mxGetPr(prhs[6]);
    double* SH = mxGetPr(prhs[7]);
    double* fc_q = mxGetPr(prhs[8]);       
        
    szA = mxGetDimensions(prhs[5]);                
    npv = szA[0];
    nc  = szA[1];
    ne  = szA[2]; 
    
    ndim = mxGetNumberOfDimensions(prhs[5]);
    if (ndim == 2) {
        ne = 1;
    }
        
    szB = mxGetDimensions(prhs[6]);                
    ndf = szB[0];
    nch = szB[1];    
    
    szC = mxGetDimensions(prhs[3]);                
    nd  = szC[2];
        
    Int ndims[5];    
    ndims[0] = nd;    
    ndims[1] = nc; 
    ndims[2] = nch;      
    ndims[3] = npv; 
    ndims[4] = ndf;    
    
    mwSize szS[3]; szS[0]=npv; szS[1]=nc; szS[2]=ne; 
    plhs[0]=mxCreateNumericArray(3,szS,mxDOUBLE_CLASS,mxREAL);         
    
    mwSize szT[5]; szT[0]=npv; szT[1]=nc; szT[2]=ndf; szT[3]=nch; szT[4]=ne; 
    plhs[1]=mxCreateNumericArray(5,szT,mxDOUBLE_CLASS,mxREAL);         
    
    double* DUDG = mxGetPr(plhs[0]);
    double* DUDG_DUH = mxGetPr(plhs[1]);
        
    Int iH, iI, iN, iO, iP, iQ, iU, iT;
    
    for (i=0; i<ne; i++) {
                
        iH = i*npv*npv*nd;
        iI = i*npv*ndf*nd;                
        iN = i*npv*nc;                
        iO = i*npv*nch;
        iP = i*npv*npv*nch*nc;        
        iQ = i*npv*nch*ndf*nch;        
        iT = i*nch*ndf;                
        iU = i*npv*nc*ndf*nch;
        
        localsolveuq(&DUDG[iN], &DUDG_DUH[iU], &Ru[iO], &BD[iP], &Rlu[iQ], &MiC[iH], &MiE[iI], 
                &UDG[iN], &UH[iT], &SH[iN], fc_q[0], ndims);                   
    }             
}
    
