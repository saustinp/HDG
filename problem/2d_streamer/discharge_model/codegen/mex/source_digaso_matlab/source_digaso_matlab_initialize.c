/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * source_digaso_matlab_initialize.c
 *
 * Code generation for function 'source_digaso_matlab_initialize'
 *
 */

/* Include files */
#include "source_digaso_matlab_initialize.h"
#include "_coder_source_digaso_matlab_mex.h"
#include "rt_nonfinite.h"
#include "source_digaso_matlab_data.h"

/* Function Declarations */
static void source_digaso_matlab_once(void);

/* Function Definitions */
static void source_digaso_matlab_once(void)
{
  mex_InitInfAndNan();
}

void source_digaso_matlab_initialize(void)
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
    source_digaso_matlab_once();
  }
}

/* End of code generation (source_digaso_matlab_initialize.c) */
