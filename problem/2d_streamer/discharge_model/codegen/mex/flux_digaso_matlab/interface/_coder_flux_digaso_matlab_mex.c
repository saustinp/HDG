/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_flux_digaso_matlab_mex.c
 *
 * Code generation for function '_coder_flux_digaso_matlab_mex'
 *
 */

/* Include files */
#include "_coder_flux_digaso_matlab_mex.h"
#include "_coder_flux_digaso_matlab_api.h"
#include "flux_digaso_matlab_data.h"
#include "flux_digaso_matlab_initialize.h"
#include "flux_digaso_matlab_terminate.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void flux_digaso_matlab_mexFunction(int32_T nlhs, mxArray *plhs[2],
                                    int32_T nrhs, const mxArray *prhs[11])
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  const mxArray *b_prhs[11];
  const mxArray *outputs[2];
  int32_T i;
  int32_T i1;
  st.tls = emlrtRootTLSGlobal;
  /* Check for proper number of arguments. */
  if (nrhs != 11) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 11, 4,
                        18, "flux_digaso_matlab");
  }
  if (nlhs > 2) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 18,
                        "flux_digaso_matlab");
  }
  /* Call the function. */
  for (i = 0; i < 11; i++) {
    b_prhs[i] = prhs[i];
  }
  flux_digaso_matlab_api(b_prhs, nlhs, outputs);
  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    i1 = 1;
  } else {
    i1 = nlhs;
  }
  emlrtReturnArrays(i1, &plhs[0], &outputs[0]);
}

void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs,
                 const mxArray *prhs[])
{
  mexAtExit(&flux_digaso_matlab_atexit);
  /* Module initialization. */
  flux_digaso_matlab_initialize();
  /* Dispatch the entry-point. */
  flux_digaso_matlab_mexFunction(nlhs, plhs, nrhs, prhs);
  /* Module termination. */
  flux_digaso_matlab_terminate();
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLSR2022a(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1,
                           NULL, "UTF-8", true);
  return emlrtRootTLSGlobal;
}

/* End of code generation (_coder_flux_digaso_matlab_mex.c) */
