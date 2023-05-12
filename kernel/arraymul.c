#include "mex.h"
#include "blas.h"

typedef ptrdiff_t  Int;

// extern "C" {
//     void dgemm(char*,char*,Int*,Int*,Int*,double*,double*,Int*,double*,Int*,double*,double*,Int*);
//     void dgemv(char*,Int*,Int*,double*,double*,Int*,double*,Int*,double*,double*,Int*);
// }

void mexFunction(int nlhs,mxArray* plhs[], int nrhs, const mxArray* prhs[])
{	
    /*#pragma omp parallel num_threads(4) */
    
	 /* Check for proper number of arguments. */
 	if(nrhs!=2) {
    	mexErrMsgTxt("two input arguments required.");
  	}        
            
    const mwSize *szA, *szB;    
    Int rowsA, colsA, rowsB, colsB, ndim3, ndimA, ndimB, i, iA, iB, iC;
    Int inc = 1;
    char *chn = "N";
    double one = 1.0, zero = 0.0;     

    /* Get dimensions of array A and array B */
    szA = mxGetDimensions(prhs[0]);
    szB = mxGetDimensions(prhs[1]);
    
    /* Get number of dimensions of array A and array B */    
    ndimA = mxGetNumberOfDimensions(prhs[0]);
    ndimB = mxGetNumberOfDimensions(prhs[1]);
    
    rowsA = szA[0];
    colsA = szA[1];
    ndim3 = szA[2];
    rowsB = szB[0];
    colsB = szB[1];


    if (ndimA!=3) {
        /* mexErrMsgTxt("A must be three-dimensional array."); */
        ndim3 = 1;
    }
    if (colsA!=rowsB) {
    	mexErrMsgTxt("Number of columns of A must be equal to numbers of rows B.");
  	}
    
    /* darray A(mxGetPr(prhs[0]),rowsA,colsA,ndim3); */     
    double* A = mxGetPr(prhs[0]);    
    double* B = mxGetPr(prhs[1]);
        
        
    if ((ndimB==2) && (colsB == ndim3)) {
        /* darray B(mxGetPr(prhs[1]),rowsB,colsB); */       
        
        mwSize szC[2]; szC[0]=rowsA; szC[1]=ndim3; 
        plhs[0]=mxCreateNumericArray(2,szC,mxDOUBLE_CLASS,mxREAL);  
        
        /* darray C(mxGetPr(plhs[0]),rowsA,ndim3); */               
        double* C = mxGetPr(plhs[0]);
        
        //#pragma omp parallel for shared(A, B, C, chn, rowsA, colsA, inc, one, zero, ndim3) private(i, iA, iB, iC) 
        for (i=0; i<ndim3; i++) {
            /* dgemv(chn, &rowsA, &colsA, &one, &A(0,0,i), &rowsA, &B(0,i), 
                    &inc, &zero, &C(0,i), &inc); */
            iA = i*rowsA*colsA;
            iB = i*rowsB;
            iC = i*rowsA;
            dgemv(chn, &rowsA, &colsA, &one, &A[iA], &rowsA, &B[iB], 
                    &inc, &zero, &C[iC], &inc);
        }        
        
    } else {
        /* darray B(mxGetPr(prhs[1]),rowsB,colsB,ndim3); */
        
        mwSize szC[3]; szC[0]=rowsA; szC[1]=colsB;  szC[2]=ndim3; 
        plhs[0]=mxCreateNumericArray(3,szC,mxDOUBLE_CLASS,mxREAL);  
        
        /* darray C(mxGetPr(plhs[0]),rowsA,colsB,ndim3); */
        double* C = mxGetPr(plhs[0]);
        
        ///#pragma omp parallel for shared(A, B, C, chn, rowsA, colsA, colsB, one, zero, ndim3) private(i, iA, iB, iC) 
        for (i=0; i<ndim3; i++) {
            /* dgemm(chn, chn, &rowsA, &colsB, &colsA, &one, &A(0,0,i), &rowsA, &B(0,0,i), 
                    &colsA, &zero, &C(0,0,i), &rowsA); */
            iA = i*rowsA*colsA;
            iB = i*rowsB*colsB;
            iC = i*rowsA*colsB;
            dgemm(chn, chn, &rowsA, &colsB, &colsA, &one, &A[iA], &rowsA, &B[iB], 
                    &colsA, &zero, &C[iC], &rowsA);
        }        
    }         
}


