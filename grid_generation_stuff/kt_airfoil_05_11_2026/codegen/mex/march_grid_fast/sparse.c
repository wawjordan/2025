/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sparse.c
 *
 * Code generation for function 'sparse'
 *
 */

/* Include files */
#include "sparse.h"
#include "eml_int_forloop_overflow_check.h"
#include "indexShapeCheck.h"
#include "march_grid_fast_data.h"
#include "march_grid_fast_emxutil.h"
#include "march_grid_fast_types.h"
#include "rt_nonfinite.h"
#include "warning.h"
#include "cs.h"
#include "makeCXSparseMatrix.h"
#include "solve_from_lu.h"
#include "solve_from_qr.h"

/* Variable Definitions */
static emlrtRSInfo eb_emlrtRSI = {
    1456,              /* lineNo */
    "sparse/mldivide", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2023a\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\sparse.m" /* pathName */
};

static emlrtRSInfo fb_emlrtRSI =
    {
        160,                         /* lineNo */
        "CXSparseAPI/iteratedSolve", /* fcnName */
        "C:\\Program "
        "Files\\MATLAB\\R2023a\\toolbox\\eml\\eml\\+coder\\+"
        "internal\\CXSparseAPI.m" /* pathName */
};

static emlrtRSInfo gb_emlrtRSI =
    {
        162,                         /* lineNo */
        "CXSparseAPI/iteratedSolve", /* fcnName */
        "C:\\Program "
        "Files\\MATLAB\\R2023a\\toolbox\\eml\\eml\\+coder\\+"
        "internal\\CXSparseAPI.m" /* pathName */
};

static emlrtRSInfo hb_emlrtRSI =
    {
        291,                      /* lineNo */
        "CXSparseAPI/iteratedLU", /* fcnName */
        "C:\\Program "
        "Files\\MATLAB\\R2023a\\toolbox\\eml\\eml\\+coder\\+"
        "internal\\CXSparseAPI.m" /* pathName */
};

static emlrtRSInfo ib_emlrtRSI =
    {
        312,                      /* lineNo */
        "CXSparseAPI/iteratedLU", /* fcnName */
        "C:\\Program "
        "Files\\MATLAB\\R2023a\\toolbox\\eml\\eml\\+coder\\+"
        "internal\\CXSparseAPI.m" /* pathName */
};

static emlrtRSInfo jb_emlrtRSI =
    {
        316,                      /* lineNo */
        "CXSparseAPI/iteratedLU", /* fcnName */
        "C:\\Program "
        "Files\\MATLAB\\R2023a\\toolbox\\eml\\eml\\+coder\\+"
        "internal\\CXSparseAPI.m" /* pathName */
};

static emlrtRSInfo kb_emlrtRSI =
    {
        457,                  /* lineNo */
        "CXSparseAPI/makeCX", /* fcnName */
        "C:\\Program "
        "Files\\MATLAB\\R2023a\\toolbox\\eml\\eml\\+coder\\+"
        "internal\\CXSparseAPI.m" /* pathName */
};

static emlrtRSInfo lb_emlrtRSI =
    {
        459,                  /* lineNo */
        "CXSparseAPI/makeCX", /* fcnName */
        "C:\\Program "
        "Files\\MATLAB\\R2023a\\toolbox\\eml\\eml\\+coder\\+"
        "internal\\CXSparseAPI.m" /* pathName */
};

static emlrtRSInfo mb_emlrtRSI = {
    395,                 /* lineNo */
    "sparse/ctranspose", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2023a\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\sparse.m" /* pathName */
};

static emlrtRSInfo nb_emlrtRSI = {
    24,                    /* lineNo */
    "sparse/locTranspose", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2023a\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\locTranspose.m" /* pathName */
};

static emlrtRSInfo ob_emlrtRSI = {
    29,                    /* lineNo */
    "sparse/locTranspose", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2023a\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\locTranspose.m" /* pathName */
};

static emlrtRSInfo pb_emlrtRSI = {
    33,                    /* lineNo */
    "sparse/locTranspose", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2023a\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\locTranspose.m" /* pathName */
};

static emlrtRSInfo qb_emlrtRSI = {
    17,                    /* lineNo */
    "sparse/locTranspose", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2023a\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\locTranspose.m" /* pathName */
};

static emlrtRSInfo tb_emlrtRSI = {
    176,             /* lineNo */
    "sparse/sparse", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2023a\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\sparse.m" /* pathName */
};

static emlrtRSInfo ub_emlrtRSI = {
    143,             /* lineNo */
    "sparse/sparse", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2023a\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\sparse.m" /* pathName */
};

static emlrtRSInfo vb_emlrtRSI = {
    142,             /* lineNo */
    "sparse/sparse", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2023a\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\sparse.m" /* pathName */
};

static emlrtRSInfo xb_emlrtRSI =
    {
        423,                      /* lineNo */
        "CXSparseAPI/iteratedQR", /* fcnName */
        "C:\\Program "
        "Files\\MATLAB\\R2023a\\toolbox\\eml\\eml\\+coder\\+"
        "internal\\CXSparseAPI.m" /* pathName */
};

static emlrtRSInfo yb_emlrtRSI =
    {
        357,                      /* lineNo */
        "CXSparseAPI/iteratedQR", /* fcnName */
        "C:\\Program "
        "Files\\MATLAB\\R2023a\\toolbox\\eml\\eml\\+coder\\+"
        "internal\\CXSparseAPI.m" /* pathName */
};

static emlrtRSInfo bc_emlrtRSI =
    {
        379,                      /* lineNo */
        "CXSparseAPI/iteratedQR", /* fcnName */
        "C:\\Program "
        "Files\\MATLAB\\R2023a\\toolbox\\eml\\eml\\+coder\\+"
        "internal\\CXSparseAPI.m" /* pathName */
};

static emlrtRSInfo cc_emlrtRSI =
    {
        380,                      /* lineNo */
        "CXSparseAPI/iteratedQR", /* fcnName */
        "C:\\Program "
        "Files\\MATLAB\\R2023a\\toolbox\\eml\\eml\\+coder\\+"
        "internal\\CXSparseAPI.m" /* pathName */
};

static emlrtMCInfo
    c_emlrtMCI =
        {
            53,        /* lineNo */
            19,        /* colNo */
            "flt2str", /* fName */
            "C:\\Program "
            "Files\\MATLAB\\R2023a\\toolbox\\eml\\eml\\+coder\\+"
            "internal\\flt2str.m" /* pName */
};

static emlrtRTEInfo o_emlrtRTEI =
    {
        139,                         /* lineNo */
        35,                          /* colNo */
        "CXSparseAPI/iteratedSolve", /* fName */
        "C:\\Program "
        "Files\\MATLAB\\R2023a\\toolbox\\eml\\eml\\+coder\\+"
        "internal\\CXSparseAPI.m" /* pName */
};

static emlrtRTEInfo p_emlrtRTEI = {
    1629,              /* lineNo */
    9,                 /* colNo */
    "assertValidSize", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2023a\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\sparse.m" /* pName */
};

static emlrtRTEInfo q_emlrtRTEI = {
    178,             /* lineNo */
    39,              /* colNo */
    "sparse/sparse", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2023a\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\sparse.m" /* pName */
};

static emlrtRTEInfo ib_emlrtRTEI = {
    1456,     /* lineNo */
    13,       /* colNo */
    "sparse", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2023a\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\sparse.m" /* pName */
};

static emlrtRTEInfo jb_emlrtRTEI =
    {
        388,           /* lineNo */
        38,            /* colNo */
        "CXSparseAPI", /* fName */
        "C:\\Program "
        "Files\\MATLAB\\R2023a\\toolbox\\eml\\eml\\+coder\\+"
        "internal\\CXSparseAPI.m" /* pName */
};

static emlrtRTEInfo kb_emlrtRTEI =
    {
        405,           /* lineNo */
        46,            /* colNo */
        "CXSparseAPI", /* fName */
        "C:\\Program "
        "Files\\MATLAB\\R2023a\\toolbox\\eml\\eml\\+coder\\+"
        "internal\\CXSparseAPI.m" /* pName */
};

static emlrtRTEInfo lb_emlrtRTEI =
    {
        399,           /* lineNo */
        46,            /* colNo */
        "CXSparseAPI", /* fName */
        "C:\\Program "
        "Files\\MATLAB\\R2023a\\toolbox\\eml\\eml\\+coder\\+"
        "internal\\CXSparseAPI.m" /* pName */
};

static emlrtRTEInfo mb_emlrtRTEI =
    {
        394,           /* lineNo */
        25,            /* colNo */
        "CXSparseAPI", /* fName */
        "C:\\Program "
        "Files\\MATLAB\\R2023a\\toolbox\\eml\\eml\\+coder\\+"
        "internal\\CXSparseAPI.m" /* pName */
};

static emlrtRTEInfo nb_emlrtRTEI =
    {
        457,           /* lineNo */
        63,            /* colNo */
        "CXSparseAPI", /* fName */
        "C:\\Program "
        "Files\\MATLAB\\R2023a\\toolbox\\eml\\eml\\+coder\\+"
        "internal\\CXSparseAPI.m" /* pName */
};

static emlrtRTEInfo ob_emlrtRTEI = {
    395,      /* lineNo */
    13,       /* colNo */
    "sparse", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2023a\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\sparse.m" /* pName */
};

static emlrtRTEInfo pb_emlrtRTEI = {
    193,      /* lineNo */
    42,       /* colNo */
    "sparse", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2023a\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\sparse.m" /* pName */
};

static emlrtRTEInfo qb_emlrtRTEI = {
    32,             /* lineNo */
    1,              /* colNo */
    "locTranspose", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2023a\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\locTranspose.m" /* pName */
};

static emlrtRSInfo
    fc_emlrtRSI =
        {
            53,        /* lineNo */
            "flt2str", /* fcnName */
            "C:\\Program "
            "Files\\MATLAB\\R2023a\\toolbox\\eml\\eml\\+coder\\+"
            "internal\\flt2str.m" /* pathName */
};

/* Function Declarations */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               char_T y[14]);

static const mxArray *b_sprintf(const emlrtStack *sp, const mxArray *m1,
                                const mxArray *m2, emlrtMCInfo *location);

static void emlrt_marshallIn(const emlrtStack *sp,
                             const mxArray *a__output_of_sprintf_,
                             const char_T *identifier, char_T y[14]);

static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                               const emlrtMsgIdentifier *msgId, char_T ret[14]);

static int32_T sparse_ctranspose(
    const emlrtStack *sp, const emxArray_real_T *this_d,
    const emxArray_int32_T *this_colidx, const emxArray_int32_T *this_rowidx,
    int32_T this_m, int32_T this_n, emxArray_real_T *y_d,
    emxArray_int32_T *y_colidx, emxArray_int32_T *y_rowidx, int32_T *y_n);

/* Function Definitions */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId, char_T y[14])
{
  g_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static const mxArray *b_sprintf(const emlrtStack *sp, const mxArray *m1,
                                const mxArray *m2, emlrtMCInfo *location)
{
  const mxArray *pArrays[2];
  const mxArray *m;
  pArrays[0] = m1;
  pArrays[1] = m2;
  return emlrtCallMATLABR2012b((emlrtConstCTX)sp, 1, &m, 2, &pArrays[0],
                               "sprintf", true, location);
}

static void emlrt_marshallIn(const emlrtStack *sp,
                             const mxArray *a__output_of_sprintf_,
                             const char_T *identifier, char_T y[14])
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  b_emlrt_marshallIn(sp, emlrtAlias(a__output_of_sprintf_), &thisId, y);
  emlrtDestroyArray(&a__output_of_sprintf_);
}

static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                               const emlrtMsgIdentifier *msgId, char_T ret[14])
{
  static const int32_T dims[2] = {1, 14};
  emlrtCheckBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "char", false, 2U,
                          (const void *)&dims[0]);
  emlrtImportCharArrayR2015b((emlrtConstCTX)sp, src, &ret[0], 14);
  emlrtDestroyArray(&src);
}

static int32_T sparse_ctranspose(
    const emlrtStack *sp, const emxArray_real_T *this_d,
    const emxArray_int32_T *this_colidx, const emxArray_int32_T *this_rowidx,
    int32_T this_m, int32_T this_n, emxArray_real_T *y_d,
    emxArray_int32_T *y_colidx, emxArray_int32_T *y_rowidx, int32_T *y_n)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack st;
  emxArray_int32_T *counts;
  const real_T *this_d_data;
  real_T *y_d_data;
  const int32_T *this_colidx_data;
  const int32_T *this_rowidx_data;
  int32_T b;
  int32_T c;
  int32_T numalloc;
  int32_T outridx;
  int32_T y_m;
  int32_T *counts_data;
  int32_T *y_colidx_data;
  int32_T *y_rowidx_data;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  this_rowidx_data = this_rowidx->data;
  this_colidx_data = this_colidx->data;
  this_d_data = this_d->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  st.site = &mb_emlrtRSI;
  b_st.site = &qb_emlrtRSI;
  c_st.site = &vb_emlrtRSI;
  if (this_n < 0) {
    emlrtErrorWithMessageIdR2018a(&c_st, &p_emlrtRTEI,
                                  "Coder:toolbox:SparseNegativeSize",
                                  "Coder:toolbox:SparseNegativeSize", 0);
  }
  if (this_n >= MAX_int32_T) {
    emlrtErrorWithMessageIdR2018a(
        &c_st, &f_emlrtRTEI, "Coder:toolbox:SparseMaxSize",
        "Coder:toolbox:SparseMaxSize", 2, 12, MAX_int32_T);
  }
  c_st.site = &ub_emlrtRSI;
  if (this_m < 0) {
    emlrtErrorWithMessageIdR2018a(&c_st, &p_emlrtRTEI,
                                  "Coder:toolbox:SparseNegativeSize",
                                  "Coder:toolbox:SparseNegativeSize", 0);
  }
  if (this_m >= MAX_int32_T) {
    emlrtErrorWithMessageIdR2018a(
        &c_st, &f_emlrtRTEI, "Coder:toolbox:SparseMaxSize",
        "Coder:toolbox:SparseMaxSize", 2, 12, MAX_int32_T);
  }
  c_st.site = &tb_emlrtRSI;
  b = this_colidx_data[this_colidx->size[0] - 1];
  if (b - 1 < 0) {
    emlrtErrorWithMessageIdR2018a(&c_st, &p_emlrtRTEI,
                                  "Coder:toolbox:SparseNegativeSize",
                                  "Coder:toolbox:SparseNegativeSize", 0);
  }
  if (b - 1 >= MAX_int32_T) {
    emlrtErrorWithMessageIdR2018a(
        &c_st, &f_emlrtRTEI, "Coder:toolbox:SparseMaxSize",
        "Coder:toolbox:SparseMaxSize", 2, 12, MAX_int32_T);
  }
  if (b - 1 < 0) {
    emlrtErrorWithMessageIdR2018a(&b_st, &q_emlrtRTEI,
                                  "Coder:toolbox:SparseNzmaxTooSmall",
                                  "Coder:toolbox:SparseNzmaxTooSmall", 0);
  }
  if (b - 1 >= 1) {
    numalloc = b - 2;
  } else {
    numalloc = 0;
  }
  outridx = y_d->size[0];
  y_d->size[0] = numalloc + 1;
  emxEnsureCapacity_real_T(&b_st, y_d, outridx, &ob_emlrtRTEI);
  y_d_data = y_d->data;
  outridx = y_colidx->size[0];
  y_colidx->size[0] = this_m + 1;
  emxEnsureCapacity_int32_T(&b_st, y_colidx, outridx, &pb_emlrtRTEI);
  y_colidx_data = y_colidx->data;
  y_colidx_data[0] = 1;
  outridx = y_rowidx->size[0];
  y_rowidx->size[0] = numalloc + 1;
  emxEnsureCapacity_int32_T(&b_st, y_rowidx, outridx, &ob_emlrtRTEI);
  y_rowidx_data = y_rowidx->data;
  for (outridx = 0; outridx <= numalloc; outridx++) {
    y_d_data[outridx] = 0.0;
    y_rowidx_data[outridx] = 0;
  }
  for (c = 0; c < this_m; c++) {
    y_colidx_data[c + 1] = 1;
  }
  outridx = y_colidx->size[0];
  for (c = 0; c <= outridx - 2; c++) {
    y_colidx_data[c] = 1;
  }
  y_colidx_data[y_colidx->size[0] - 1] = 1;
  emxInit_int32_T(&st, &counts, &qb_emlrtRTEI);
  if ((this_m != 0) && (this_n != 0)) {
    numalloc = y_colidx->size[0];
    outridx = y_colidx->size[0];
    y_colidx->size[0] = numalloc;
    emxEnsureCapacity_int32_T(&st, y_colidx, outridx, &ob_emlrtRTEI);
    y_colidx_data = y_colidx->data;
    for (outridx = 0; outridx < numalloc; outridx++) {
      y_colidx_data[outridx] = 0;
    }
    b_st.site = &nb_emlrtRSI;
    if (b - 1 > 2147483646) {
      c_st.site = &w_emlrtRSI;
      check_forloop_overflow_error(&c_st);
    }
    for (numalloc = 0; numalloc <= b - 2; numalloc++) {
      y_colidx_data[this_rowidx_data[numalloc]]++;
    }
    y_colidx_data[0] = 1;
    b = this_m + 1;
    b_st.site = &ob_emlrtRSI;
    if (this_m + 1 > 2147483646) {
      c_st.site = &w_emlrtRSI;
      check_forloop_overflow_error(&c_st);
    }
    for (numalloc = 2; numalloc <= b; numalloc++) {
      y_colidx_data[numalloc - 1] += y_colidx_data[numalloc - 2];
    }
    outridx = counts->size[0];
    counts->size[0] = this_m;
    emxEnsureCapacity_int32_T(&st, counts, outridx, &qb_emlrtRTEI);
    counts_data = counts->data;
    for (outridx = 0; outridx < this_m; outridx++) {
      counts_data[outridx] = 0;
    }
    b_st.site = &pb_emlrtRSI;
    for (c = 0; c < this_n; c++) {
      for (numalloc = this_colidx_data[c] - 1;
           numalloc + 1 < this_colidx_data[c + 1]; numalloc++) {
        b = counts_data[this_rowidx_data[numalloc] - 1];
        outridx = (b + y_colidx_data[this_rowidx_data[numalloc] - 1]) - 1;
        y_d_data[outridx] = this_d_data[numalloc];
        y_rowidx_data[outridx] = c + 1;
        counts_data[this_rowidx_data[numalloc] - 1] = b + 1;
      }
    }
  }
  emxFree_int32_T(&st, &counts);
  y_m = this_n;
  *y_n = this_m;
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
  return y_m;
}

void sparse_mldivide(const emlrtStack *sp, const emxArray_real_T *A_d,
                     const emxArray_int32_T *A_colidx,
                     const emxArray_int32_T *A_rowidx, int32_T A_m, int32_T A_n,
                     const emxArray_real_T *b, emxArray_real_T *y)
{
  static const int32_T iv[2] = {1, 6};
  static const char_T rfmt[6] = {'%', '1', '4', '.', '6', 'e'};
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack st;
  emxArray_int32_T *in_colidx;
  emxArray_int32_T *in_rowidx;
  emxArray_real_T *outBuff;
  const mxArray *b_y;
  const mxArray *c_y;
  const mxArray *m;
  const real_T *A_d_data;
  const real_T *b_data;
  real_T *outBuff_data;
  real_T *y_data;
  const int32_T *A_colidx_data;
  const int32_T *A_rowidx_data;
  int32_T in_n;
  int32_T *in_colidx_data;
  int32_T *in_rowidx_data;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  e_st.prev = &d_st;
  e_st.tls = d_st.tls;
  f_st.prev = &e_st;
  f_st.tls = e_st.tls;
  b_data = b->data;
  A_rowidx_data = A_rowidx->data;
  A_colidx_data = A_colidx->data;
  A_d_data = A_d->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  st.site = &eb_emlrtRSI;
  emxInit_real_T(&st, &outBuff, 1, &mb_emlrtRTEI);
  emxInit_int32_T(&st, &in_colidx, &nb_emlrtRTEI);
  emxInit_int32_T(&st, &in_rowidx, &nb_emlrtRTEI);
  if ((A_m == 0) || (A_n == 0) || (b->size[0] == 0)) {
    in_n = y->size[0];
    y->size[0] = A_n;
    emxEnsureCapacity_real_T(&st, y, in_n, &ib_emlrtRTEI);
    y_data = y->data;
    for (in_n = 0; in_n < A_n; in_n++) {
      y_data[in_n] = 0.0;
    }
  } else {
    if (A_m != b->size[0]) {
      emlrtErrorWithMessageIdR2018a(&st, &o_emlrtRTEI, "MATLAB:dimagree",
                                    "MATLAB:dimagree", 0);
    }
    if (b->size[0] == A_n) {
      cs_di *b_cxA;
      cs_din *b_N;
      cs_dis *b_S;
      int32_T rank;
      b_st.site = &fb_emlrtRSI;
      c_st.site = &hb_emlrtRSI;
      if (A_m < A_n) {
        d_st.site = &kb_emlrtRSI;
        e_st.site = &kb_emlrtRSI;
        rank = sparse_ctranspose(&e_st, A_d, A_colidx, A_rowidx, A_m, A_n,
                                 outBuff, in_colidx, in_rowidx, &in_n);
        in_rowidx_data = in_rowidx->data;
        in_colidx_data = in_colidx->data;
        outBuff_data = outBuff->data;
        b_cxA = makeCXSparseMatrix(in_colidx_data[in_colidx->size[0] - 1] - 1,
                                   in_n, rank, &in_colidx_data[0],
                                   &in_rowidx_data[0], &outBuff_data[0]);
      } else {
        d_st.site = &lb_emlrtRSI;
        b_cxA = makeCXSparseMatrix(A_colidx_data[A_colidx->size[0] - 1] - 1,
                                   A_n, A_m, &A_colidx_data[0],
                                   &A_rowidx_data[0], &A_d_data[0]);
      }
      b_S = cs_di_sqr(2, b_cxA, 0);
      b_N = cs_di_lu(b_cxA, b_S, 1);
      cs_di_spfree(b_cxA);
      if (b_N == NULL) {
        cs_di *c_cxA;
        cs_din *c_N;
        cs_dis *c_S;
        real_T tol;
        c_st.site = &ib_emlrtRSI;
        warning(&c_st);
        cs_di_sfree(b_S);
        cs_di_nfree(b_N);
        c_st.site = &jb_emlrtRSI;
        d_st.site = &yb_emlrtRSI;
        if (A_m < A_n) {
          e_st.site = &kb_emlrtRSI;
          f_st.site = &kb_emlrtRSI;
          rank = sparse_ctranspose(&f_st, A_d, A_colidx, A_rowidx, A_m, A_n,
                                   outBuff, in_colidx, in_rowidx, &in_n);
          in_rowidx_data = in_rowidx->data;
          in_colidx_data = in_colidx->data;
          outBuff_data = outBuff->data;
          c_cxA = makeCXSparseMatrix(in_colidx_data[in_colidx->size[0] - 1] - 1,
                                     in_n, rank, &in_colidx_data[0],
                                     &in_rowidx_data[0], &outBuff_data[0]);
        } else {
          e_st.site = &lb_emlrtRSI;
          c_cxA = makeCXSparseMatrix(A_colidx_data[A_colidx->size[0] - 1] - 1,
                                     A_n, A_m, &A_colidx_data[0],
                                     &A_rowidx_data[0], &A_d_data[0]);
        }
        c_S = cs_di_sqr(2, c_cxA, 1);
        c_N = cs_di_qr(c_cxA, c_S);
        cs_di_spfree(c_cxA);
        qr_rank_di(c_N, &tol);
        in_n = y->size[0];
        y->size[0] = A_n;
        emxEnsureCapacity_real_T(&c_st, y, in_n, &jb_emlrtRTEI);
        y_data = y->data;
        if (b->size[0] < A_n) {
          in_n = outBuff->size[0];
          outBuff->size[0] = A_n;
          emxEnsureCapacity_real_T(&c_st, outBuff, in_n, &lb_emlrtRTEI);
          outBuff_data = outBuff->data;
        } else {
          in_n = outBuff->size[0];
          outBuff->size[0] = b->size[0];
          emxEnsureCapacity_real_T(&c_st, outBuff, in_n, &kb_emlrtRTEI);
          outBuff_data = outBuff->data;
        }
        rank = b->size[0];
        for (in_n = 0; in_n < rank; in_n++) {
          outBuff_data[in_n] = b_data[in_n];
        }
        int32_T iv1[2];
        solve_from_qr_di(c_N, c_S, (double *)&outBuff_data[0], b->size[0], A_n);
        iv1[0] = 1;
        iv1[1] = A_n;
        d_st.site = &xb_emlrtRSI;
        indexShapeCheck(&d_st, outBuff->size[0], iv1);
        for (in_n = 0; in_n < A_n; in_n++) {
          y_data[in_n] = outBuff_data[in_n];
        }
        cs_di_sfree(c_S);
        cs_di_nfree(c_N);
      } else {
        in_n = y->size[0];
        y->size[0] = b->size[0];
        emxEnsureCapacity_real_T(&b_st, y, in_n, &ib_emlrtRTEI);
        y_data = y->data;
        rank = b->size[0];
        for (in_n = 0; in_n < rank; in_n++) {
          y_data[in_n] = b_data[in_n];
        }
        solve_from_lu_di(b_N, b_S, (double *)&y_data[0], b->size[0]);
        cs_di_sfree(b_S);
        cs_di_nfree(b_N);
      }
    } else {
      cs_di *cxA;
      cs_din *N;
      cs_dis *S;
      real_T tol;
      int32_T iv1[2];
      int32_T rank;
      b_st.site = &gb_emlrtRSI;
      c_st.site = &yb_emlrtRSI;
      if (A_m < A_n) {
        d_st.site = &kb_emlrtRSI;
        e_st.site = &kb_emlrtRSI;
        rank = sparse_ctranspose(&e_st, A_d, A_colidx, A_rowidx, A_m, A_n,
                                 outBuff, in_colidx, in_rowidx, &in_n);
        in_rowidx_data = in_rowidx->data;
        in_colidx_data = in_colidx->data;
        outBuff_data = outBuff->data;
        cxA = makeCXSparseMatrix(in_colidx_data[in_colidx->size[0] - 1] - 1,
                                 in_n, rank, &in_colidx_data[0],
                                 &in_rowidx_data[0], &outBuff_data[0]);
      } else {
        d_st.site = &lb_emlrtRSI;
        cxA = makeCXSparseMatrix(A_colidx_data[A_colidx->size[0] - 1] - 1, A_n,
                                 A_m, &A_colidx_data[0], &A_rowidx_data[0],
                                 &A_d_data[0]);
      }
      S = cs_di_sqr(2, cxA, 1);
      N = cs_di_qr(cxA, S);
      cs_di_spfree(cxA);
      rank = qr_rank_di(N, &tol);
      if (A_m > A_n) {
        in_n = A_n;
      } else {
        in_n = A_m;
      }
      if (rank < in_n) {
        char_T str[14];
        c_st.site = &cc_emlrtRSI;
        b_y = NULL;
        m = emlrtCreateCharArray(2, &iv[0]);
        emlrtInitCharArrayR2013a(&c_st, 6, m, &rfmt[0]);
        emlrtAssign(&b_y, m);
        c_y = NULL;
        m = emlrtCreateDoubleScalar(tol);
        emlrtAssign(&c_y, m);
        d_st.site = &fc_emlrtRSI;
        emlrt_marshallIn(&d_st, b_sprintf(&d_st, b_y, c_y, &c_emlrtMCI),
                         "<output of sprintf>", str);
        c_st.site = &bc_emlrtRSI;
        b_warning(&c_st, rank, str);
      }
      in_n = y->size[0];
      y->size[0] = A_n;
      emxEnsureCapacity_real_T(&b_st, y, in_n, &jb_emlrtRTEI);
      y_data = y->data;
      if (b->size[0] < A_n) {
        in_n = outBuff->size[0];
        outBuff->size[0] = A_n;
        emxEnsureCapacity_real_T(&b_st, outBuff, in_n, &lb_emlrtRTEI);
        outBuff_data = outBuff->data;
      } else {
        in_n = outBuff->size[0];
        outBuff->size[0] = b->size[0];
        emxEnsureCapacity_real_T(&b_st, outBuff, in_n, &kb_emlrtRTEI);
        outBuff_data = outBuff->data;
      }
      rank = b->size[0];
      for (in_n = 0; in_n < rank; in_n++) {
        outBuff_data[in_n] = b_data[in_n];
      }
      solve_from_qr_di(N, S, (double *)&outBuff_data[0], b->size[0], A_n);
      if (A_n < 1) {
        rank = 0;
      } else {
        rank = A_n;
      }
      iv1[0] = 1;
      iv1[1] = rank;
      c_st.site = &xb_emlrtRSI;
      indexShapeCheck(&c_st, outBuff->size[0], iv1);
      for (in_n = 0; in_n < rank; in_n++) {
        y_data[in_n] = outBuff_data[in_n];
      }
      cs_di_sfree(S);
      cs_di_nfree(N);
    }
  }
  emxFree_int32_T(&st, &in_rowidx);
  emxFree_int32_T(&st, &in_colidx);
  emxFree_real_T(&st, &outBuff);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

/* End of code generation (sparse.c) */
