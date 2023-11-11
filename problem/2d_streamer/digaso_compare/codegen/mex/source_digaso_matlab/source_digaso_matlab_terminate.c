/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * source_digaso_matlab_terminate.c
 *
 * Code generation for function 'source_digaso_matlab_terminate'
 *
 */

/* Include files */
#include "source_digaso_matlab_terminate.h"
#include "_coder_source_digaso_matlab_mex.h"
#include "rt_nonfinite.h"
#include "source_digaso_matlab_data.h"

/* Function Definitions */
void source_digaso_matlab_atexit(void)
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

void source_digaso_matlab_terminate(void)
{
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (source_digaso_matlab_terminate.c) */
