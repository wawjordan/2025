/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_march_grid_fast_api.c
 *
 * Code generation for function '_coder_march_grid_fast_api'
 *
 */

/* Include files */
#include "_coder_march_grid_fast_api.h"
#include "march_grid_fast.h"
#include "march_grid_fast_data.h"
#include "march_grid_fast_emxutil.h"
#include "march_grid_fast_types.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRTEInfo rb_emlrtRTEI = {
    1,                            /* lineNo */
    1,                            /* colNo */
    "_coder_march_grid_fast_api", /* fName */
    ""                            /* pName */
};

/* Function Declarations */
static real_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *imax,
                                 const char_T *identifier);

static real_T d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                 const emlrtMsgIdentifier *parentId);

static void e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *x,
                               const char_T *identifier, emxArray_real_T *y);

static const mxArray *emlrt_marshallOut(const emxArray_real_T *u);

static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               emxArray_real_T *y);

static real_T h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                 const emlrtMsgIdentifier *msgId);

static void i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               emxArray_real_T *ret);

/* Function Definitions */
static real_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *imax,
                                 const char_T *identifier)
{
  emlrtMsgIdentifier thisId;
  real_T y;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = d_emlrt_marshallIn(sp, emlrtAlias(imax), &thisId);
  emlrtDestroyArray(&imax);
  return y;
}

static real_T d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                 const emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = h_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static void e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *x,
                               const char_T *identifier, emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  f_emlrt_marshallIn(sp, emlrtAlias(x), &thisId, y);
  emlrtDestroyArray(&x);
}

static const mxArray *emlrt_marshallOut(const emxArray_real_T *u)
{
  static const int32_T i = 0;
  const mxArray *m;
  const mxArray *y;
  const real_T *u_data;
  u_data = u->data;
  y = NULL;
  m = emlrtCreateNumericArray(1, (const void *)&i, mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m, (void *)&u_data[0]);
  emlrtSetDimensions((mxArray *)m, &u->size[0], 1);
  emlrtAssign(&y, m);
  return y;
}

static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               emxArray_real_T *y)
{
  i_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static real_T h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
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

static void i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               emxArray_real_T *ret)
{
  static const int32_T dims = 8192;
  int32_T i;
  int32_T i1;
  boolean_T b = true;
  emlrtCheckVsBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 1U,
                            (const void *)&dims, &b, &i);
  ret->allocatedSize = i;
  i1 = ret->size[0];
  ret->size[0] = i;
  emxEnsureCapacity_real_T(sp, ret, i1, (emlrtRTEInfo *)NULL);
  ret->data = (real_T *)emlrtMxGetData(src);
  ret->canFreeData = false;
  emlrtDestroyArray(&src);
}

void march_grid_fast_api(const mxArray *const prhs[10], int32_T nlhs,
                         const mxArray *plhs[2])
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  emxArray_real_T *alpham;
  emxArray_real_T *mu;
  emxArray_real_T *muim;
  emxArray_real_T *rj;
  emxArray_real_T *rjm1;
  emxArray_real_T *scale;
  emxArray_real_T *x;
  emxArray_real_T *x_update;
  emxArray_real_T *y;
  emxArray_real_T *y_update;
  real_T imax;
  real_T jmax;
  st.tls = emlrtRootTLSGlobal;
  emlrtHeapReferenceStackEnterFcnR2012b(&st);
  /* Marshall function inputs */
  imax = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[0]), "imax");
  jmax = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "jmax");
  emxInit_real_T(&st, &x, 1, &rb_emlrtRTEI);
  x->canFreeData = false;
  e_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "x", x);
  emxInit_real_T(&st, &y, 1, &rb_emlrtRTEI);
  y->canFreeData = false;
  e_emlrt_marshallIn(&st, emlrtAlias(prhs[3]), "y", y);
  emxInit_real_T(&st, &mu, 1, &rb_emlrtRTEI);
  mu->canFreeData = false;
  e_emlrt_marshallIn(&st, emlrtAlias(prhs[4]), "mu", mu);
  emxInit_real_T(&st, &muim, 1, &rb_emlrtRTEI);
  muim->canFreeData = false;
  e_emlrt_marshallIn(&st, emlrtAlias(prhs[5]), "muim", muim);
  emxInit_real_T(&st, &alpham, 1, &rb_emlrtRTEI);
  alpham->canFreeData = false;
  e_emlrt_marshallIn(&st, emlrtAlias(prhs[6]), "alpham", alpham);
  emxInit_real_T(&st, &scale, 1, &rb_emlrtRTEI);
  scale->canFreeData = false;
  e_emlrt_marshallIn(&st, emlrtAlias(prhs[7]), "scale", scale);
  emxInit_real_T(&st, &rj, 1, &rb_emlrtRTEI);
  rj->canFreeData = false;
  e_emlrt_marshallIn(&st, emlrtAlias(prhs[8]), "rj", rj);
  emxInit_real_T(&st, &rjm1, 1, &rb_emlrtRTEI);
  rjm1->canFreeData = false;
  e_emlrt_marshallIn(&st, emlrtAlias(prhs[9]), "rjm1", rjm1);
  /* Invoke the target function */
  emxInit_real_T(&st, &x_update, 1, &rb_emlrtRTEI);
  emxInit_real_T(&st, &y_update, 1, &rb_emlrtRTEI);
  march_grid_fast(&st, imax, jmax, x, y, mu, muim, alpham, scale, rj, rjm1,
                  x_update, y_update);
  emxFree_real_T(&st, &rjm1);
  emxFree_real_T(&st, &rj);
  emxFree_real_T(&st, &scale);
  emxFree_real_T(&st, &alpham);
  emxFree_real_T(&st, &muim);
  emxFree_real_T(&st, &mu);
  emxFree_real_T(&st, &y);
  emxFree_real_T(&st, &x);
  /* Marshall function outputs */
  x_update->canFreeData = false;
  plhs[0] = emlrt_marshallOut(x_update);
  emxFree_real_T(&st, &x_update);
  if (nlhs > 1) {
    y_update->canFreeData = false;
    plhs[1] = emlrt_marshallOut(y_update);
  }
  emxFree_real_T(&st, &y_update);
  emlrtHeapReferenceStackLeaveFcnR2012b(&st);
}

/* End of code generation (_coder_march_grid_fast_api.c) */
