/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * fhat_digaso_matlab_initialize.c
 *
 * Code generation for function 'fhat_digaso_matlab_initialize'
 *
 */

/* Include files */
#include "fhat_digaso_matlab_initialize.h"
#include "_coder_fhat_digaso_matlab_mex.h"
#include "fhat_digaso_matlab_data.h"
#include "rt_nonfinite.h"

/* Function Declarations */
static void fhat_digaso_matlab_once(void);

/* Function Definitions */
static void fhat_digaso_matlab_once(void)
{
  mex_InitInfAndNan();
}

void fhat_digaso_matlab_initialize(void)
{
  static const volatile char_T *emlrtBreakCheckR2012bFlagVar = NULL;
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtBreakCheckR2012bFlagVar = emlrtGetBreakCheckFlagAddressR2022b(&st);
  emlrtClearAllocCountR2012b(&st, false, 0U, NULL);
  emlrtEnterRtStackR2012b(&st);
  if (emlrtFirstTimeR2012b(emlrtRootTLSGlobal)) {
    fhat_digaso_matlab_once();
  }
}

/* End of code generation (fhat_digaso_matlab_initialize.c) */
