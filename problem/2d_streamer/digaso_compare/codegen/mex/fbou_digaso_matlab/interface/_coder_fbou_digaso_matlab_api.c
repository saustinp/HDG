/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_fbou_digaso_matlab_api.c
 *
 * Code generation for function '_coder_fbou_digaso_matlab_api'
 *
 */

/* Include files */
#include "_coder_fbou_digaso_matlab_api.h"
#include "fbou_digaso_matlab.h"
#include "fbou_digaso_matlab_data.h"
#include "rt_nonfinite.h"

/* Function Declarations */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[234];

static void b_emlrt_marshallOut(const real_T u[2106], const mxArray *y);

static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *fh_u,
                                   const char_T *identifier))[2106];

static void c_emlrt_marshallOut(const real_T u[702], const mxArray *y);

static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[2106];

static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *fh_uh,
                                   const char_T *identifier))[702];

static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *fh,
                                 const char_T *identifier))[234];

static void emlrt_marshallOut(const real_T u[234], const mxArray *y);

static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[702];

static real_T (*g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *pg,
                                   const char_T *identifier))[156];

static real_T (*h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[156];

static real_T (*i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *ui,
                                   const char_T *identifier))[78];

static real_T (*j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[78];

static real_T (*k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *param,
                                   const char_T *identifier))[10];

static real_T (*l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[10];

static real_T m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *b_time,
                                 const char_T *identifier);

static real_T n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                 const emlrtMsgIdentifier *parentId);

static real_T (*o_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[234];

static real_T (*p_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[2106];

static real_T (*q_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[702];

static real_T (*r_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[156];

static real_T (*s_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[78];

static real_T (*t_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[10];

static real_T u_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                 const emlrtMsgIdentifier *msgId);

/* Function Definitions */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[234]
{
  real_T(*y)[234];
  y = o_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static void b_emlrt_marshallOut(const real_T u[2106], const mxArray *y)
{
  static const int32_T iv[2] = {1, 2106};
  emlrtMxSetData((mxArray *)y, (void *)&u[0]);
  emlrtSetDimensions((mxArray *)y, &iv[0], 2);
}

static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *fh_u,
                                   const char_T *identifier))[2106]
{
  emlrtMsgIdentifier thisId;
  real_T(*y)[2106];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = d_emlrt_marshallIn(sp, emlrtAlias(fh_u), &thisId);
  emlrtDestroyArray(&fh_u);
  return y;
}

static void c_emlrt_marshallOut(const real_T u[702], const mxArray *y)
{
  static const int32_T iv[2] = {1, 702};
  emlrtMxSetData((mxArray *)y, (void *)&u[0]);
  emlrtSetDimensions((mxArray *)y, &iv[0], 2);
}

static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[2106]
{
  real_T(*y)[2106];
  y = p_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *fh_uh,
                                   const char_T *identifier))[702]
{
  emlrtMsgIdentifier thisId;
  real_T(*y)[702];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = f_emlrt_marshallIn(sp, emlrtAlias(fh_uh), &thisId);
  emlrtDestroyArray(&fh_uh);
  return y;
}

static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *fh,
                                 const char_T *identifier))[234]
{
  emlrtMsgIdentifier thisId;
  real_T(*y)[234];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(fh), &thisId);
  emlrtDestroyArray(&fh);
  return y;
}

static void emlrt_marshallOut(const real_T u[234], const mxArray *y)
{
  static const int32_T iv[2] = {1, 234};
  emlrtMxSetData((mxArray *)y, (void *)&u[0]);
  emlrtSetDimensions((mxArray *)y, &iv[0], 2);
}

static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[702]
{
  real_T(*y)[702];
  y = q_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *pg,
                                   const char_T *identifier))[156]
{
  emlrtMsgIdentifier thisId;
  real_T(*y)[156];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = h_emlrt_marshallIn(sp, emlrtAlias(pg), &thisId);
  emlrtDestroyArray(&pg);
  return y;
}

static real_T (*h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[156]
{
  real_T(*y)[156];
  y = r_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *ui,
                                   const char_T *identifier))[78]
{
  emlrtMsgIdentifier thisId;
  real_T(*y)[78];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = j_emlrt_marshallIn(sp, emlrtAlias(ui), &thisId);
  emlrtDestroyArray(&ui);
  return y;
}

static real_T (*j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[78]
{
  real_T(*y)[78];
  y = s_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *param,
                                   const char_T *identifier))[10]
{
  emlrtMsgIdentifier thisId;
  real_T(*y)[10];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = l_emlrt_marshallIn(sp, emlrtAlias(param), &thisId);
  emlrtDestroyArray(&param);
  return y;
}

static real_T (*l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[10]
{
  real_T(*y)[10];
  y = t_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *b_time,
                                 const char_T *identifier)
{
  emlrtMsgIdentifier thisId;
  real_T y;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = n_emlrt_marshallIn(sp, emlrtAlias(b_time), &thisId);
  emlrtDestroyArray(&b_time);
  return y;
}

static real_T n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                 const emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = u_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*o_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[234]
{
  static const int32_T dims[2] = {1, 234};
  real_T(*ret)[234];
  int32_T iv[2];
  boolean_T bv[2] = {false, false};
  emlrtCheckVsBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 2U,
                            (const void *)&dims[0], &bv[0], &iv[0]);
  ret = (real_T(*)[234])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T (*p_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[2106]
{
  static const int32_T dims[2] = {1, 2106};
  real_T(*ret)[2106];
  int32_T iv[2];
  boolean_T bv[2] = {false, false};
  emlrtCheckVsBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 2U,
                            (const void *)&dims[0], &bv[0], &iv[0]);
  ret = (real_T(*)[2106])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T (*q_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[702]
{
  static const int32_T dims[2] = {1, 702};
  real_T(*ret)[702];
  int32_T iv[2];
  boolean_T bv[2] = {false, false};
  emlrtCheckVsBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 2U,
                            (const void *)&dims[0], &bv[0], &iv[0]);
  ret = (real_T(*)[702])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T (*r_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[156]
{
  static const int32_T dims[2] = {1, 156};
  real_T(*ret)[156];
  int32_T iv[2];
  boolean_T bv[2] = {false, false};
  emlrtCheckVsBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 2U,
                            (const void *)&dims[0], &bv[0], &iv[0]);
  ret = (real_T(*)[156])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T (*s_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[78]
{
  static const int32_T dims[2] = {1, 78};
  real_T(*ret)[78];
  int32_T iv[2];
  boolean_T bv[2] = {false, false};
  emlrtCheckVsBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 2U,
                            (const void *)&dims[0], &bv[0], &iv[0]);
  ret = (real_T(*)[78])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T (*t_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
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

static real_T u_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
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

void fbou_digaso_matlab_api(const mxArray *const prhs[16], int32_T nlhs,
                            const mxArray *plhs[3])
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
  const mxArray *prhs_copy_idx_5;
  const mxArray *prhs_copy_idx_6;
  const mxArray *prhs_copy_idx_7;
  const mxArray *prhs_copy_idx_8;
  real_T(*fh_u)[2106];
  real_T(*fh_uh)[702];
  real_T(*udg)[702];
  real_T(*fh)[234];
  real_T(*uhg)[234];
  real_T(*nl)[156];
  real_T(*pg)[156];
  real_T(*ui)[78];
  real_T(*param)[10];
  real_T b_time;
  real_T ib;
  real_T nc;
  real_T ncd;
  real_T nch;
  real_T nd;
  real_T ng;
  st.tls = emlrtRootTLSGlobal;
  prhs_copy_idx_0 = emlrtProtectR2012b(prhs[0], 0, true, -1);
  prhs_copy_idx_1 = emlrtProtectR2012b(prhs[1], 1, true, -1);
  prhs_copy_idx_2 = emlrtProtectR2012b(prhs[2], 2, true, -1);
  prhs_copy_idx_3 = emlrtProtectR2012b(prhs[3], 3, false, -1);
  prhs_copy_idx_4 = emlrtProtectR2012b(prhs[4], 4, false, -1);
  prhs_copy_idx_5 = emlrtProtectR2012b(prhs[5], 5, false, -1);
  prhs_copy_idx_6 = emlrtProtectR2012b(prhs[6], 6, false, -1);
  prhs_copy_idx_7 = emlrtProtectR2012b(prhs[7], 7, false, -1);
  prhs_copy_idx_8 = emlrtProtectR2012b(prhs[8], 8, false, -1);
  /* Marshall function inputs */
  fh = emlrt_marshallIn(&st, emlrtAlias(prhs_copy_idx_0), "fh");
  fh_u = c_emlrt_marshallIn(&st, emlrtAlias(prhs_copy_idx_1), "fh_u");
  fh_uh = e_emlrt_marshallIn(&st, emlrtAlias(prhs_copy_idx_2), "fh_uh");
  pg = g_emlrt_marshallIn(&st, emlrtAlias(prhs_copy_idx_3), "pg");
  udg = e_emlrt_marshallIn(&st, emlrtAlias(prhs_copy_idx_4), "udg");
  uhg = emlrt_marshallIn(&st, emlrtAlias(prhs_copy_idx_5), "uhg");
  nl = g_emlrt_marshallIn(&st, emlrtAlias(prhs_copy_idx_6), "nl");
  ui = i_emlrt_marshallIn(&st, emlrtAlias(prhs_copy_idx_7), "ui");
  param = k_emlrt_marshallIn(&st, emlrtAlias(prhs_copy_idx_8), "param");
  b_time = m_emlrt_marshallIn(&st, emlrtAliasP(prhs[9]), "time");
  ib = m_emlrt_marshallIn(&st, emlrtAliasP(prhs[10]), "ib");
  ng = m_emlrt_marshallIn(&st, emlrtAliasP(prhs[11]), "ng");
  nc = m_emlrt_marshallIn(&st, emlrtAliasP(prhs[12]), "nc");
  nch = m_emlrt_marshallIn(&st, emlrtAliasP(prhs[13]), "nch");
  nd = m_emlrt_marshallIn(&st, emlrtAliasP(prhs[14]), "nd");
  ncd = m_emlrt_marshallIn(&st, emlrtAliasP(prhs[15]), "ncd");
  /* Invoke the target function */
  fbou_digaso_matlab(&st, *fh, *fh_u, *fh_uh, *pg, *udg, *uhg, *nl, *ui, *param,
                     b_time, ib, ng, nc, nch, nd, ncd);
  /* Marshall function outputs */
  emlrt_marshallOut(*fh, prhs_copy_idx_0);
  plhs[0] = prhs_copy_idx_0;
  if (nlhs > 1) {
    b_emlrt_marshallOut(*fh_u, prhs_copy_idx_1);
    plhs[1] = prhs_copy_idx_1;
  }
  if (nlhs > 2) {
    c_emlrt_marshallOut(*fh_uh, prhs_copy_idx_2);
    plhs[2] = prhs_copy_idx_2;
  }
}

/* End of code generation (_coder_fbou_digaso_matlab_api.c) */
