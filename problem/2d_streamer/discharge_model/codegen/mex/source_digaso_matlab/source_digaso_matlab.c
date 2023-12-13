/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * source_digaso_matlab.c
 *
 * Code generation for function 'source_digaso_matlab'
 *
 */

/* Include files */
#include "source_digaso_matlab.h"
#include "rt_nonfinite.h"
#include "mwmathutil.h"
#include "source.h"

/* Function Definitions */
void source_digaso_matlab(const emlrtStack *sp, real_T s[90000],
                          real_T s_udg[810000], real_T pg[60000],
                          real_T udg[270000], real_T param[10], real_T b_time,
                          real_T ng, real_T nc, real_T ncu, real_T nd,
                          real_T ncd)
{
  real_T d;
  int32_T i;
  int32_T i1;
  int32_T i2;
  int32_T i3;
  int32_T i4;
  int32_T i5;
  (void)sp;
  /*  Note that the arrays are pre-flattened before passing into the */
  /*  function */
  /*  Include the header for the C code */
  /*  Evaluate the C function - the source function returns void */
  d = muDoubleScalarRound(b_time);
  if (d < 2.147483648E+9) {
    if (d >= -2.147483648E+9) {
      i = (int32_T)d;
    } else {
      i = MIN_int32_T;
    }
  } else if (d >= 2.147483648E+9) {
    i = MAX_int32_T;
  } else {
    i = 0;
  }
  d = muDoubleScalarRound(ng);
  if (d < 2.147483648E+9) {
    if (d >= -2.147483648E+9) {
      i1 = (int32_T)d;
    } else {
      i1 = MIN_int32_T;
    }
  } else if (d >= 2.147483648E+9) {
    i1 = MAX_int32_T;
  } else {
    i1 = 0;
  }
  d = muDoubleScalarRound(nc);
  if (d < 2.147483648E+9) {
    if (d >= -2.147483648E+9) {
      i2 = (int32_T)d;
    } else {
      i2 = MIN_int32_T;
    }
  } else if (d >= 2.147483648E+9) {
    i2 = MAX_int32_T;
  } else {
    i2 = 0;
  }
  d = muDoubleScalarRound(ncu);
  if (d < 2.147483648E+9) {
    if (d >= -2.147483648E+9) {
      i3 = (int32_T)d;
    } else {
      i3 = MIN_int32_T;
    }
  } else if (d >= 2.147483648E+9) {
    i3 = MAX_int32_T;
  } else {
    i3 = 0;
  }
  d = muDoubleScalarRound(nd);
  if (d < 2.147483648E+9) {
    if (d >= -2.147483648E+9) {
      i4 = (int32_T)d;
    } else {
      i4 = MIN_int32_T;
    }
  } else if (d >= 2.147483648E+9) {
    i4 = MAX_int32_T;
  } else {
    i4 = 0;
  }
  d = muDoubleScalarRound(ncd);
  if (d < 2.147483648E+9) {
    if (d >= -2.147483648E+9) {
      i5 = (int32_T)d;
    } else {
      i5 = MIN_int32_T;
    }
  } else if (d >= 2.147483648E+9) {
    i5 = MAX_int32_T;
  } else {
    i5 = 0;
  }
  source(&s[0], &s_udg[0], &pg[0], &udg[0], &param[0], i, i1, i2, i3, i4, i5);
}

/* End of code generation (source_digaso_matlab.c) */
