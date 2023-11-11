/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_source_digaso_matlab_api.c
 *
 * Code generation for function '_coder_source_digaso_matlab_api'
 *
 */

/* Include files */
#include "_coder_source_digaso_matlab_api.h"
#include "rt_nonfinite.h"
#include "source_digaso_matlab.h"
#include "source_digaso_matlab_data.h"

/* Function Declarations */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[90000];

static void b_emlrt_marshallOut(const real_T u[810000], const mxArray *y);

static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *s_udg,
                                   const char_T *identifier))[810000];

static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[810000];

static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *pg,
                                   const char_T *identifier))[60000];

static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *s,
                                 const char_T *identifier))[90000];

static void emlrt_marshallOut(const real_T u[90000], const mxArray *y);

static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[60000];

static real_T (*g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *udg,
                                   const char_T *identifier))[270000];

static real_T (*h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[270000];

static real_T (*i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *param,
                                   const char_T *identifier))[10];

static real_T (*j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[10];

static real_T k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *b_time,
                                 const char_T *identifier);

static real_T l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                 const emlrtMsgIdentifier *parentId);

static real_T (*m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[90000];

static real_T (*n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[810000];

static real_T (*o_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[60000];

static real_T (*p_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[270000];

static real_T (*q_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[10];

static real_T r_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                 const emlrtMsgIdentifier *msgId);

/* Function Definitions */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[90000]
{
  real_T(*y)[90000];
  y = m_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static void b_emlrt_marshallOut(const real_T u[810000], const mxArray *y)
{
  static const int32_T iv[2] = {1, 810000};
  emlrtMxSetData((mxArray *)y, (void *)&u[0]);
  emlrtSetDimensions((mxArray *)y, &iv[0], 2);
}

static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *s_udg,
                                   const char_T *identifier))[810000]
{
  emlrtMsgIdentifier thisId;
  real_T(*y)[810000];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = d_emlrt_marshallIn(sp, emlrtAlias(s_udg), &thisId);
  emlrtDestroyArray(&s_udg);
  return y;
}

static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[810000]
{
  real_T(*y)[810000];
  y = n_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *pg,
                                   const char_T *identifier))[60000]
{
  emlrtMsgIdentifier thisId;
  real_T(*y)[60000];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = f_emlrt_marshallIn(sp, emlrtAlias(pg), &thisId);
  emlrtDestroyArray(&pg);
  return y;
}

static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *s,
                                 const char_T *identifier))[90000]
{
  emlrtMsgIdentifier thisId;
  real_T(*y)[90000];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(s), &thisId);
  emlrtDestroyArray(&s);
  return y;
}

static void emlrt_marshallOut(const real_T u[90000], const mxArray *y)
{
  static const int32_T iv[2] = {1, 90000};
  emlrtMxSetData((mxArray *)y, (void *)&u[0]);
  emlrtSetDimensions((mxArray *)y, &iv[0], 2);
}

static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[60000]
{
  real_T(*y)[60000];
  y = o_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *udg,
                                   const char_T *identifier))[270000]
{
  emlrtMsgIdentifier thisId;
  real_T(*y)[270000];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = h_emlrt_marshallIn(sp, emlrtAlias(udg), &thisId);
  emlrtDestroyArray(&udg);
  return y;
}

static real_T (*h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[270000]
{
  real_T(*y)[270000];
  y = p_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *param,
                                   const char_T *identifier))[10]
{
  emlrtMsgIdentifier thisId;
  real_T(*y)[10];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = j_emlrt_marshallIn(sp, emlrtAlias(param), &thisId);
  emlrtDestroyArray(&param);
  return y;
}

static real_T (*j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[10]
{
  real_T(*y)[10];
  y = q_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *b_time,
                                 const char_T *identifier)
{
  emlrtMsgIdentifier thisId;
  real_T y;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = l_emlrt_marshallIn(sp, emlrtAlias(b_time), &thisId);
  emlrtDestroyArray(&b_time);
  return y;
}

static real_T l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                 const emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = r_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[90000]
{
  static const int32_T dims[2] = {1, 90000};
  real_T(*ret)[90000];
  int32_T iv[2];
  boolean_T bv[2] = {false, false};
  emlrtCheckVsBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 2U,
                            (const void *)&dims[0], &bv[0], &iv[0]);
  ret = (real_T(*)[90000])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T (*n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[810000]
{
  static const int32_T dims[2] = {1, 810000};
  real_T(*ret)[810000];
  int32_T iv[2];
  boolean_T bv[2] = {false, false};
  emlrtCheckVsBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 2U,
                            (const void *)&dims[0], &bv[0], &iv[0]);
  ret = (real_T(*)[810000])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T (*o_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[60000]
{
  static const int32_T dims[2] = {1, 60000};
  real_T(*ret)[60000];
  int32_T iv[2];
  boolean_T bv[2] = {false, false};
  emlrtCheckVsBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 2U,
                            (const void *)&dims[0], &bv[0], &iv[0]);
  ret = (real_T(*)[60000])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T (*p_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[270000]
{
  static const int32_T dims[2] = {1, 270000};
  real_T(*ret)[270000];
  int32_T iv[2];
  boolean_T bv[2] = {false, false};
  emlrtCheckVsBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 2U,
                            (const void *)&dims[0], &bv[0], &iv[0]);
  ret = (real_T(*)[270000])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T (*q_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[10]
{
  static const int32_T dims[2] = {1, 10};
  real_T(*ret)[10];
  int32_T iv[2];
  boolean_T bv[2] = {false, false};
  emlrtCheckVsBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 2U,
                            (const void *)&dims[0], &bv[0], &iv[0]);
  ret = (real_T(*)[10])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T r_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                 const emlrtMsgIdentifier *msgId)
{
  static const int32_T dims = 0;
  real_T ret;
  emlrtCheckBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 0U,
                          (const void *)&dims);
  ret = *(real_T *)emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

void source_digaso_matlab_api(const mxArray *const prhs[11], int32_T nlhs,
                              const mxArray *plhs[2])
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  const mxArray *prhs_copy_idx_0;
  const mxArray *prhs_copy_idx_1;
  const mxArray *prhs_copy_idx_2;
  const mxArray *prhs_copy_idx_3;
  const mxArray *prhs_copy_idx_4;
  real_T(*s_udg)[810000];
  real_T(*udg)[270000];
  real_T(*s)[90000];
  real_T(*pg)[60000];
  real_T(*param)[10];
  real_T b_time;
  real_T nc;
  real_T ncd;
  real_T ncu;
  real_T nd;
  real_T ng;
  st.tls = emlrtRootTLSGlobal;
  prhs_copy_idx_0 = emlrtProtectR2012b(prhs[0], 0, true, -1);
  prhs_copy_idx_1 = emlrtProtectR2012b(prhs[1], 1, true, -1);
  prhs_copy_idx_2 = emlrtProtectR2012b(prhs[2], 2, false, -1);
  prhs_copy_idx_3 = emlrtProtectR2012b(prhs[3], 3, false, -1);
  prhs_copy_idx_4 = emlrtProtectR2012b(prhs[4], 4, false, -1);
  /* Marshall function inputs */
  s = emlrt_marshallIn(&st, emlrtAlias(prhs_copy_idx_0), "s");
  s_udg = c_emlrt_marshallIn(&st, emlrtAlias(prhs_copy_idx_1), "s_udg");
  pg = e_emlrt_marshallIn(&st, emlrtAlias(prhs_copy_idx_2), "pg");
  udg = g_emlrt_marshallIn(&st, emlrtAlias(prhs_copy_idx_3), "udg");
  param = i_emlrt_marshallIn(&st, emlrtAlias(prhs_copy_idx_4), "param");
  b_time = k_emlrt_marshallIn(&st, emlrtAliasP(prhs[5]), "time");
  ng = k_emlrt_marshallIn(&st, emlrtAliasP(prhs[6]), "ng");
  nc = k_emlrt_marshallIn(&st, emlrtAliasP(prhs[7]), "nc");
  ncu = k_emlrt_marshallIn(&st, emlrtAliasP(prhs[8]), "ncu");
  nd = k_emlrt_marshallIn(&st, emlrtAliasP(prhs[9]), "nd");
  ncd = k_emlrt_marshallIn(&st, emlrtAliasP(prhs[10]), "ncd");
  /* Invoke the target function */
  source_digaso_matlab(&st, *s, *s_udg, *pg, *udg, *param, b_time, ng, nc, ncu,
                       nd, ncd);
  /* Marshall function outputs */
  emlrt_marshallOut(*s, prhs_copy_idx_0);
  plhs[0] = prhs_copy_idx_0;
  if (nlhs > 1) {
    b_emlrt_marshallOut(*s_udg, prhs_copy_idx_1);
    plhs[1] = prhs_copy_idx_1;
  }
}

/* End of code generation (_coder_source_digaso_matlab_api.c) */
