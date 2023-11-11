/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * fhat_digaso_matlab_terminate.c
 *
 * Code generation for function 'fhat_digaso_matlab_terminate'
 *
 */

/* Include files */
#include "fhat_digaso_matlab_terminate.h"
#include "_coder_fhat_digaso_matlab_mex.h"
#include "fhat_digaso_matlab_data.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void fhat_digaso_matlab_atexit(void)
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

void fhat_digaso_matlab_terminate(void)
{
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (fhat_digaso_matlab_terminate.c) */
