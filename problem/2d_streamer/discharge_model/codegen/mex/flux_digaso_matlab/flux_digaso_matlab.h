/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * flux_digaso_matlab.h
 *
 * Code generation for function 'flux_digaso_matlab'
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
void flux_digaso_matlab(const emlrtStack *sp, real_T f[180000],
                        real_T f_udg[1620000], real_T pg[60000],
                        real_T udg[270000], real_T param[10], real_T b_time,
                        real_T ng, real_T nc, real_T ncu, real_T nd,
                        real_T ncd);

/* End of code generation (flux_digaso_matlab.h) */
