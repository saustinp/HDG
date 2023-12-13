/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * fbou_digaso_matlab_terminate.c
 *
 * Code generation for function 'fbou_digaso_matlab_terminate'
 *
 */

/* Include files */
#include "fbou_digaso_matlab_terminate.h"
#include "_coder_fbou_digaso_matlab_mex.h"
#include "fbou_digaso_matlab_data.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void fbou_digaso_matlab_atexit(void)
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

void fbou_digaso_matlab_terminate(void)
{
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (fbou_digaso_matlab_terminate.c) */
