/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * flux_digaso_matlab_terminate.c
 *
 * Code generation for function 'flux_digaso_matlab_terminate'
 *
 */

/* Include files */
#include "flux_digaso_matlab_terminate.h"
#include "_coder_flux_digaso_matlab_mex.h"
#include "flux_digaso_matlab_data.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void flux_digaso_matlab_atexit(void)
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

void flux_digaso_matlab_terminate(void)
{
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (flux_digaso_matlab_terminate.c) */
