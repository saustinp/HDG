/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_mathOpsIntegrated_mex.c
 *
 * Code generation for function '_coder_mathOpsIntegrated_mex'
 *
 */

/* Include files */
#include "_coder_mathOpsIntegrated_mex.h"
#include "_coder_mathOpsIntegrated_api.h"
#include "mathOpsIntegrated_data.h"
#include "mathOpsIntegrated_initialize.h"
#include "mathOpsIntegrated_terminate.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void mathOpsIntegrated_mexFunction(int32_T nlhs, mxArray *plhs[2], int32_T nrhs,
                                   const mxArray *prhs[2])
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  const mxArray *outputs[2];
  int32_T i;
  st.tls = emlrtRootTLSGlobal;
  /* Check for proper number of arguments. */
  if (nrhs != 2) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 2, 4,
                        17, "mathOpsIntegrated");
  }
  if (nlhs > 2) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 17,
                        "mathOpsIntegrated");
  }
  /* Call the function. */
  mathOpsIntegrated_api(prhs, nlhs, outputs);
  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    i = 1;
  } else {
    i = nlhs;
  }
  emlrtReturnArrays(i, &plhs[0], &outputs[0]);
}

void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs,
                 const mxArray *prhs[])
{
  mexAtExit(&mathOpsIntegrated_atexit);
  /* Module initialization. */
  mathOpsIntegrated_initialize();
  /* Dispatch the entry-point. */
  mathOpsIntegrated_mexFunction(nlhs, plhs, nrhs, prhs);
  /* Module termination. */
  mathOpsIntegrated_terminate();
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLSR2022a(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1,
                           NULL, "UTF-8", true);
  return emlrtRootTLSGlobal;
}

/* End of code generation (_coder_mathOpsIntegrated_mex.c) */
