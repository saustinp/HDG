/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * fbou_digaso_matlab.h
 *
 * Code generation for function 'fbou_digaso_matlab'
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
void fbou_digaso_matlab(const emlrtStack *sp, real_T fh[198], real_T fh_u[1782],
                        real_T fh_uh[594], real_T pg[132], real_T udg[594],
                        real_T uhg[198], real_T nl[132], real_T ui[66],
                        real_T param[10], real_T b_time, real_T ib, real_T ng,
                        real_T nc, real_T nch, real_T nd, real_T ncd);

/* End of code generation (fbou_digaso_matlab.h) */
