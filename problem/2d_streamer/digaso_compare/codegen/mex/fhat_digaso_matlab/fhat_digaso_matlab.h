/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * fhat_digaso_matlab.h
 *
 * Code generation for function 'fhat_digaso_matlab'
 *
 */

#pragma once

/* Include files */
#include "rtwtypes.h"
#include "emlrt.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
void fhat_digaso_matlab(const emlrtStack *sp, real_T fh[135000],
                        real_T fh_udg[1215000], real_T fh_uh[405000],
                        real_T pg[90000], real_T udg[405000], real_T uh[135000],
                        real_T nl[90000], real_T param[10], real_T b_time,
                        real_T ng, real_T nc, real_T nch, real_T nd,
                        real_T ncd);

/* End of code generation (fhat_digaso_matlab.h) */
