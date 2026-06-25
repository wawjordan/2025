/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * march_grid_fast.c
 *
 * Code generation for function 'march_grid_fast'
 *
 */

/* Include files */
#include "march_grid_fast.h"
#include "assertValidSizeArg.h"
#include "eml_int_forloop_overflow_check.h"
#include "march_grid_fast_data.h"
#include "march_grid_fast_emxutil.h"
#include "march_grid_fast_types.h"
#include "rt_nonfinite.h"
#include "sparse.h"
#include "mwmathutil.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = {
    20,                /* lineNo */
    "march_grid_fast", /* fcnName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m" /* pathName */
};

static emlrtRSInfo b_emlrtRSI = {
    23,                /* lineNo */
    "march_grid_fast", /* fcnName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m" /* pathName */
};

static emlrtRSInfo c_emlrtRSI = {
    26,                /* lineNo */
    "march_grid_fast", /* fcnName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m" /* pathName */
};

static emlrtRSInfo d_emlrtRSI = {
    27,                /* lineNo */
    "march_grid_fast", /* fcnName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m" /* pathName */
};

static emlrtRSInfo e_emlrtRSI = {
    30,                /* lineNo */
    "march_grid_fast", /* fcnName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m" /* pathName */
};

static emlrtRSInfo f_emlrtRSI = {
    32,                /* lineNo */
    "march_grid_fast", /* fcnName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m" /* pathName */
};

static emlrtRSInfo g_emlrtRSI = {
    33,                /* lineNo */
    "march_grid_fast", /* fcnName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m" /* pathName */
};

static emlrtRSInfo h_emlrtRSI = {
    34,                /* lineNo */
    "march_grid_fast", /* fcnName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m" /* pathName */
};

static emlrtRSInfo i_emlrtRSI = {
    44,            /* lineNo */
    "calc_volume", /* fcnName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m" /* pathName */
};

static emlrtRSInfo j_emlrtRSI = {
    49,            /* lineNo */
    "calc_volume", /* fcnName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m" /* pathName */
};

static emlrtRSInfo k_emlrtRSI = {
    44,       /* lineNo */
    "mpower", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2023a\\toolbox\\eml\\lib\\matlab\\matfun\\mpower.m" /* pathName
                                                                          */
};

static emlrtRSInfo l_emlrtRSI =
    {
        71,      /* lineNo */
        "power", /* fcnName */
        "C:\\Program "
        "Files\\MATLAB\\R2023a\\toolbox\\eml\\lib\\matlab\\ops\\power.m" /* pathName
                                                                          */
};

static emlrtRSInfo m_emlrtRSI = {
    75,             /* lineNo */
    "grid_metrics", /* fcnName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m" /* pathName */
};

static emlrtRSInfo o_emlrtRSI = {
    91,             /* lineNo */
    "grid_metrics", /* fcnName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m" /* pathName */
};

static emlrtRSInfo q_emlrtRSI = {
    151,          /* lineNo */
    "create_lhs", /* fcnName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m" /* pathName */
};

static emlrtRSInfo r_emlrtRSI = {
    159,          /* lineNo */
    "create_lhs", /* fcnName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m" /* pathName */
};

static emlrtRSInfo s_emlrtRSI = {
    223,          /* lineNo */
    "create_lhs", /* fcnName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m" /* pathName */
};

static emlrtRSInfo t_emlrtRSI =
    {
        50,    /* lineNo */
        "eye", /* fcnName */
        "C:\\Program "
        "Files\\MATLAB\\R2023a\\toolbox\\eml\\lib\\matlab\\elmat\\eye.m" /* pathName
                                                                          */
};

static emlrtRSInfo u_emlrtRSI =
    {
        96,    /* lineNo */
        "eye", /* fcnName */
        "C:\\Program "
        "Files\\MATLAB\\R2023a\\toolbox\\eml\\lib\\matlab\\elmat\\eye.m" /* pathName
                                                                          */
};

static emlrtRSInfo v_emlrtRSI = {
    21,                           /* lineNo */
    "checkAndSaturateExpandSize", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2023a\\toolbox\\eml\\eml\\+coder\\+"
    "internal\\checkAndSaturateExpandSize.m" /* pathName */
};

static emlrtRSInfo x_emlrtRSI = {
    13,       /* lineNo */
    "sparse", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2023a\\toolbox\\eml\\lib\\matlab\\sparfun\\sparse.m" /* pathName
                                                                           */
};

static emlrtRSInfo y_emlrtRSI = {
    49,              /* lineNo */
    "sparse/sparse", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2023a\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\sparse.m" /* pathName */
};

static emlrtRSInfo ab_emlrtRSI = {
    50,              /* lineNo */
    "sparse/sparse", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2023a\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\sparse.m" /* pathName */
};

static emlrtRSInfo bb_emlrtRSI = {
    65,              /* lineNo */
    "sparse/sparse", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2023a\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\sparse.m" /* pathName */
};

static emlrtRSInfo cb_emlrtRSI = {
    66,              /* lineNo */
    "sparse/sparse", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2023a\\toolbox\\shared\\coder\\coder\\+coder\\+internal\\@"
    "sparse\\sparse.m" /* pathName */
};

static emlrtRSInfo dc_emlrtRSI = {
    40,                  /* lineNo */
    "reshapeSizeChecks", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2023a\\toolbox\\eml\\eml\\+coder\\+"
    "internal\\reshapeSizeChecks.m" /* pathName */
};

static emlrtBCInfo emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    130,          /* lineNo */
    7,            /* colNo */
    "rhs",        /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo b_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    129,          /* lineNo */
    7,            /* colNo */
    "rhs",        /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo c_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    128,          /* lineNo */
    42,           /* colNo */
    "alpha",      /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo d_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    128,          /* lineNo */
    28,           /* colNo */
    "volume",     /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo e_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    128,          /* lineNo */
    18,           /* colNo */
    "alpha",      /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo f_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    117,          /* lineNo */
    95,           /* colNo */
    "y",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo g_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    117,          /* lineNo */
    83,           /* colNo */
    "y",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo h_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    117,          /* lineNo */
    72,           /* colNo */
    "y",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo i_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    117,          /* lineNo */
    60,           /* colNo */
    "y",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo j_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    117,          /* lineNo */
    48,           /* colNo */
    "y",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo k_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    117,          /* lineNo */
    36,           /* colNo */
    "y",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo l_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    117,          /* lineNo */
    27,           /* colNo */
    "mu",         /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo m_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    120,          /* lineNo */
    83,           /* colNo */
    "y",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo n_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    120,          /* lineNo */
    72,           /* colNo */
    "y",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo o_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    120,          /* lineNo */
    60,           /* colNo */
    "y",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo p_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    120,          /* lineNo */
    48,           /* colNo */
    "y",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo q_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    120,          /* lineNo */
    36,           /* colNo */
    "y",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo r_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    120,          /* lineNo */
    27,           /* colNo */
    "mu",         /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo s_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    142,          /* lineNo */
    5,            /* colNo */
    "rhs",        /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo t_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    116,          /* lineNo */
    95,           /* colNo */
    "x",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo u_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    116,          /* lineNo */
    83,           /* colNo */
    "x",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo v_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    116,          /* lineNo */
    72,           /* colNo */
    "x",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo w_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    116,          /* lineNo */
    60,           /* colNo */
    "x",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo x_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    116,          /* lineNo */
    48,           /* colNo */
    "x",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo y_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    116,          /* lineNo */
    36,           /* colNo */
    "x",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo ab_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    116,          /* lineNo */
    27,           /* colNo */
    "mu",         /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo bb_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    119,          /* lineNo */
    83,           /* colNo */
    "x",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo cb_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    119,          /* lineNo */
    72,           /* colNo */
    "x",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo db_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    119,          /* lineNo */
    60,           /* colNo */
    "x",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo eb_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    119,          /* lineNo */
    48,           /* colNo */
    "x",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo fb_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    119,          /* lineNo */
    36,           /* colNo */
    "x",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo gb_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    119,          /* lineNo */
    27,           /* colNo */
    "mu",         /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo hb_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    141,          /* lineNo */
    5,            /* colNo */
    "rhs",        /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo ib_emlrtBCI = {
    -1,            /* iFirst */
    -1,            /* iLast */
    50,            /* lineNo */
    10,            /* colNo */
    "volume",      /* aName */
    "calc_volume", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo jb_emlrtBCI = {
    -1,            /* iFirst */
    -1,            /* iLast */
    50,            /* lineNo */
    45,            /* colNo */
    "scale",       /* aName */
    "calc_volume", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo kb_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    111,          /* lineNo */
    95,           /* colNo */
    "y",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo lb_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    111,          /* lineNo */
    84,           /* colNo */
    "y",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo mb_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    111,          /* lineNo */
    72,           /* colNo */
    "y",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo nb_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    111,          /* lineNo */
    60,           /* colNo */
    "y",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo ob_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    111,          /* lineNo */
    48,           /* colNo */
    "y",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo pb_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    111,          /* lineNo */
    36,           /* colNo */
    "y",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo qb_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    111,          /* lineNo */
    27,           /* colNo */
    "mu",         /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo rb_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    110,          /* lineNo */
    95,           /* colNo */
    "x",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo sb_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    110,          /* lineNo */
    84,           /* colNo */
    "x",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo tb_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    110,          /* lineNo */
    72,           /* colNo */
    "x",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo ub_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    110,          /* lineNo */
    60,           /* colNo */
    "x",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo vb_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    110,          /* lineNo */
    48,           /* colNo */
    "x",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo wb_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    110,          /* lineNo */
    36,           /* colNo */
    "x",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo xb_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    110,          /* lineNo */
    27,           /* colNo */
    "mu",         /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo yb_emlrtBCI = {
    -1,            /* iFirst */
    -1,            /* iLast */
    49,            /* lineNo */
    57,            /* colNo */
    "y",           /* aName */
    "calc_volume", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo ac_emlrtBCI = {
    -1,            /* iFirst */
    -1,            /* iLast */
    49,            /* lineNo */
    48,            /* colNo */
    "y",           /* aName */
    "calc_volume", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtDCInfo emlrtDCI = {
    26,                /* lineNo */
    1,                 /* colNo */
    "march_grid_fast", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    1                                                  /* checkKind */
};

static emlrtBCInfo bc_emlrtBCI = {
    -1,            /* iFirst */
    -1,            /* iLast */
    49,            /* lineNo */
    28,            /* colNo */
    "x",           /* aName */
    "calc_volume", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo cc_emlrtBCI = {
    -1,            /* iFirst */
    -1,            /* iLast */
    49,            /* lineNo */
    19,            /* colNo */
    "x",           /* aName */
    "calc_volume", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo dc_emlrtBCI = {
    -1,            /* iFirst */
    -1,            /* iLast */
    44,            /* lineNo */
    65,            /* colNo */
    "y",           /* aName */
    "calc_volume", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo ec_emlrtBCI = {
    -1,            /* iFirst */
    -1,            /* iLast */
    44,            /* lineNo */
    58,            /* colNo */
    "y",           /* aName */
    "calc_volume", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo fc_emlrtBCI = {
    -1,            /* iFirst */
    -1,            /* iLast */
    57,            /* lineNo */
    8,             /* colNo */
    "volume",      /* aName */
    "calc_volume", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo gc_emlrtBCI = {
    -1,            /* iFirst */
    -1,            /* iLast */
    57,            /* lineNo */
    20,            /* colNo */
    "volume",      /* aName */
    "calc_volume", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo hc_emlrtBCI = {
    -1,            /* iFirst */
    -1,            /* iLast */
    44,            /* lineNo */
    47,            /* colNo */
    "x",           /* aName */
    "calc_volume", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo ic_emlrtBCI = {
    -1,            /* iFirst */
    -1,            /* iLast */
    44,            /* lineNo */
    40,            /* colNo */
    "x",           /* aName */
    "calc_volume", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtDCInfo b_emlrtDCI = {
    20,                /* lineNo */
    1,                 /* colNo */
    "march_grid_fast", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    4                                                  /* checkKind */
};

static emlrtDCInfo c_emlrtDCI = {
    20,                /* lineNo */
    1,                 /* colNo */
    "march_grid_fast", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    1                                                  /* checkKind */
};

static emlrtECInfo emlrtECI = {
    1,                 /* nDims */
    33,                /* lineNo */
    12,                /* colNo */
    "march_grid_fast", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m" /* pName */
};

static emlrtECInfo b_emlrtECI = {
    1,                 /* nDims */
    34,                /* lineNo */
    12,                /* colNo */
    "march_grid_fast", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m" /* pName */
};

static emlrtRTEInfo emlrtRTEI = {
    74,                  /* lineNo */
    13,                  /* colNo */
    "reshapeSizeChecks", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2023a\\toolbox\\eml\\eml\\+coder\\+"
    "internal\\reshapeSizeChecks.m" /* pName */
};

static emlrtRTEInfo b_emlrtRTEI = {
    81,                  /* lineNo */
    23,                  /* colNo */
    "reshapeSizeChecks", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2023a\\toolbox\\eml\\eml\\+coder\\+"
    "internal\\reshapeSizeChecks.m" /* pName */
};

static emlrtRTEInfo c_emlrtRTEI = {
    43,            /* lineNo */
    9,             /* colNo */
    "calc_volume", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m" /* pName */
};

static emlrtRTEInfo d_emlrtRTEI = {
    47,            /* lineNo */
    9,             /* colNo */
    "calc_volume", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m" /* pName */
};

static emlrtBCInfo jc_emlrtBCI = {
    -1,            /* iFirst */
    -1,            /* iLast */
    48,            /* lineNo */
    13,            /* colNo */
    "rj",          /* aName */
    "calc_volume", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo kc_emlrtBCI = {
    -1,            /* iFirst */
    -1,            /* iLast */
    48,            /* lineNo */
    25,            /* colNo */
    "rjm1",        /* aName */
    "calc_volume", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo lc_emlrtBCI = {
    -1,            /* iFirst */
    -1,            /* iLast */
    48,            /* lineNo */
    35,            /* colNo */
    "rj",          /* aName */
    "calc_volume", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo mc_emlrtBCI = {
    -1,            /* iFirst */
    -1,            /* iLast */
    48,            /* lineNo */
    46,            /* colNo */
    "rjm1",        /* aName */
    "calc_volume", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo nc_emlrtBCI = {
    -1,            /* iFirst */
    -1,            /* iLast */
    54,            /* lineNo */
    20,            /* colNo */
    "volume",      /* aName */
    "calc_volume", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo oc_emlrtBCI = {
    -1,            /* iFirst */
    -1,            /* iLast */
    54,            /* lineNo */
    8,             /* colNo */
    "volume",      /* aName */
    "calc_volume", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo pc_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    107,          /* lineNo */
    27,           /* colNo */
    "mu",         /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo qc_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    107,          /* lineNo */
    36,           /* colNo */
    "x",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo rc_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    107,          /* lineNo */
    48,           /* colNo */
    "x",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo sc_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    107,          /* lineNo */
    60,           /* colNo */
    "x",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo tc_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    107,          /* lineNo */
    72,           /* colNo */
    "x",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo uc_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    107,          /* lineNo */
    84,           /* colNo */
    "x",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo vc_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    107,          /* lineNo */
    95,           /* colNo */
    "x",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo wc_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    108,          /* lineNo */
    27,           /* colNo */
    "mu",         /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo xc_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    108,          /* lineNo */
    36,           /* colNo */
    "y",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo yc_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    108,          /* lineNo */
    48,           /* colNo */
    "y",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo ad_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    108,          /* lineNo */
    60,           /* colNo */
    "y",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo bd_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    108,          /* lineNo */
    72,           /* colNo */
    "y",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo cd_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    108,          /* lineNo */
    84,           /* colNo */
    "y",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo dd_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    108,          /* lineNo */
    95,           /* colNo */
    "y",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo ed_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    113,          /* lineNo */
    27,           /* colNo */
    "mu",         /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo fd_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    113,          /* lineNo */
    36,           /* colNo */
    "x",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo gd_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    113,          /* lineNo */
    48,           /* colNo */
    "x",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo hd_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    113,          /* lineNo */
    60,           /* colNo */
    "x",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo id_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    113,          /* lineNo */
    72,           /* colNo */
    "x",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo jd_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    113,          /* lineNo */
    83,           /* colNo */
    "x",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo kd_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    113,          /* lineNo */
    95,           /* colNo */
    "x",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo ld_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    114,          /* lineNo */
    27,           /* colNo */
    "mu",         /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo md_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    114,          /* lineNo */
    36,           /* colNo */
    "y",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo nd_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    114,          /* lineNo */
    48,           /* colNo */
    "y",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo od_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    114,          /* lineNo */
    60,           /* colNo */
    "y",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo pd_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    114,          /* lineNo */
    72,           /* colNo */
    "y",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo qd_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    114,          /* lineNo */
    83,           /* colNo */
    "y",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo rd_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    114,          /* lineNo */
    95,           /* colNo */
    "y",          /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo sd_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    122,          /* lineNo */
    15,           /* colNo */
    "xxi",        /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo td_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    123,          /* lineNo */
    15,           /* colNo */
    "yxi",        /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo ud_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    124,          /* lineNo */
    16,           /* colNo */
    "xeta",       /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo vd_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    125,          /* lineNo */
    16,           /* colNo */
    "yeta",       /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo wd_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    137,          /* lineNo */
    5,            /* colNo */
    "rhs",        /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo xd_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    138,          /* lineNo */
    5,            /* colNo */
    "rhs",        /* aName */
    "create_rhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtRTEInfo e_emlrtRTEI = {
    69,             /* lineNo */
    9,              /* colNo */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m" /* pName */
};

static emlrtBCInfo yd_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    75,             /* lineNo */
    18,             /* colNo */
    "xxi",          /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo ae_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    75,             /* lineNo */
    29,             /* colNo */
    "yxi",          /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo be_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    82,             /* lineNo */
    21,             /* colNo */
    "x",            /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo ce_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    82,             /* lineNo */
    30,             /* colNo */
    "x",            /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo de_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    82,             /* lineNo */
    39,             /* colNo */
    "x",            /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo ee_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    82,             /* lineNo */
    5,              /* colNo */
    "xxi",          /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo fe_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    83,             /* lineNo */
    21,             /* colNo */
    "y",            /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo ge_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    83,             /* lineNo */
    30,             /* colNo */
    "y",            /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo he_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    83,             /* lineNo */
    39,             /* colNo */
    "y",            /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo ie_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    83,             /* lineNo */
    5,              /* colNo */
    "yxi",          /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo je_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    84,             /* lineNo */
    14,             /* colNo */
    "xxi",          /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo ke_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    84,             /* lineNo */
    25,             /* colNo */
    "yxi",          /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo le_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    85,             /* lineNo */
    16,             /* colNo */
    "yxi",          /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo me_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    85,             /* lineNo */
    26,             /* colNo */
    "volume",       /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo ne_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    86,             /* lineNo */
    16,             /* colNo */
    "xxi",          /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo oe_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    86,             /* lineNo */
    26,             /* colNo */
    "volume",       /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo pe_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    91,             /* lineNo */
    14,             /* colNo */
    "xxi",          /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo qe_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    91,             /* lineNo */
    25,             /* colNo */
    "yxi",          /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo re_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    85,             /* lineNo */
    6,              /* colNo */
    "xeta",         /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo se_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    86,             /* lineNo */
    6,              /* colNo */
    "yeta",         /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtDCInfo d_emlrtDCI = {
    64,             /* lineNo */
    1,              /* colNo */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    1                                                  /* checkKind */
};

static emlrtDCInfo e_emlrtDCI = {
    64,             /* lineNo */
    1,              /* colNo */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    4                                                  /* checkKind */
};

static emlrtDCInfo f_emlrtDCI = {
    65,             /* lineNo */
    1,              /* colNo */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    1                                                  /* checkKind */
};

static emlrtDCInfo g_emlrtDCI = {
    66,             /* lineNo */
    1,              /* colNo */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    1                                                  /* checkKind */
};

static emlrtDCInfo h_emlrtDCI = {
    67,             /* lineNo */
    1,              /* colNo */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    1                                                  /* checkKind */
};

static emlrtBCInfo te_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    71,             /* lineNo */
    22,             /* colNo */
    "x",            /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo ue_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    71,             /* lineNo */
    31,             /* colNo */
    "x",            /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo ve_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    71,             /* lineNo */
    9,              /* colNo */
    "xxi",          /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo we_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    72,             /* lineNo */
    22,             /* colNo */
    "y",            /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo xe_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    72,             /* lineNo */
    31,             /* colNo */
    "y",            /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo ye_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    72,             /* lineNo */
    9,              /* colNo */
    "yxi",          /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo af_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    77,             /* lineNo */
    20,             /* colNo */
    "yxi",          /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo bf_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    77,             /* lineNo */
    30,             /* colNo */
    "volume",       /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo cf_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    77,             /* lineNo */
    10,             /* colNo */
    "xeta",         /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo df_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    78,             /* lineNo */
    20,             /* colNo */
    "xxi",          /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo ef_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    78,             /* lineNo */
    30,             /* colNo */
    "volume",       /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo ff_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    78,             /* lineNo */
    10,             /* colNo */
    "yeta",         /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo gf_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    89,             /* lineNo */
    21,             /* colNo */
    "x",            /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo hf_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    89,             /* lineNo */
    30,             /* colNo */
    "x",            /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo if_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    89,             /* lineNo */
    39,             /* colNo */
    "x",            /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo jf_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    89,             /* lineNo */
    5,              /* colNo */
    "xxi",          /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo kf_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    90,             /* lineNo */
    21,             /* colNo */
    "y",            /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo lf_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    90,             /* lineNo */
    30,             /* colNo */
    "y",            /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo mf_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    90,             /* lineNo */
    39,             /* colNo */
    "y",            /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo nf_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    90,             /* lineNo */
    5,              /* colNo */
    "yxi",          /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo of_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    92,             /* lineNo */
    26,             /* colNo */
    "volume",       /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo pf_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    92,             /* lineNo */
    6,              /* colNo */
    "xeta",         /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo qf_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    93,             /* lineNo */
    26,             /* colNo */
    "volume",       /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo rf_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    93,             /* lineNo */
    6,              /* colNo */
    "yeta",         /* aName */
    "grid_metrics", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtRTEInfo g_emlrtRTEI = {
    156,          /* lineNo */
    9,            /* colNo */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m" /* pName */
};

static emlrtBCInfo sf_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    157,          /* lineNo */
    14,           /* colNo */
    "xxi",        /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo tf_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    158,          /* lineNo */
    14,           /* colNo */
    "yxi",        /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo uf_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    166,          /* lineNo */
    43,           /* colNo */
    "alpha",      /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtRTEInfo h_emlrtRTEI = {
    176,          /* lineNo */
    9,            /* colNo */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m" /* pName */
};

static emlrtBCInfo vf_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    183,          /* lineNo */
    32,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo wf_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    183,          /* lineNo */
    38,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo xf_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    183,          /* lineNo */
    54,           /* colNo */
    "muim",       /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo yf_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    183,          /* lineNo */
    9,            /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo ag_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    183,          /* lineNo */
    15,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo bg_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    184,          /* lineNo */
    32,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo cg_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    184,          /* lineNo */
    38,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo dg_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    184,          /* lineNo */
    54,           /* colNo */
    "muim",       /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo eg_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    184,          /* lineNo */
    9,            /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo fg_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    184,          /* lineNo */
    15,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtRTEInfo i_emlrtRTEI = {
    185,          /* lineNo */
    9,            /* colNo */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m" /* pName */
};

static emlrtBCInfo gg_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    192,          /* lineNo */
    32,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo hg_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    192,          /* lineNo */
    38,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo ig_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    192,          /* lineNo */
    54,           /* colNo */
    "muim",       /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo jg_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    192,          /* lineNo */
    9,            /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo kg_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    192,          /* lineNo */
    15,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo lg_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    193,          /* lineNo */
    32,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo mg_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    193,          /* lineNo */
    38,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo ng_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    193,          /* lineNo */
    54,           /* colNo */
    "muim",       /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo og_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    193,          /* lineNo */
    9,            /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo pg_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    193,          /* lineNo */
    15,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtRTEInfo j_emlrtRTEI = {
    194,          /* lineNo */
    9,            /* colNo */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m" /* pName */
};

static emlrtRTEInfo k_emlrtRTEI = {
    203,          /* lineNo */
    9,            /* colNo */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m" /* pName */
};

static emlrtRTEInfo l_emlrtRTEI = {
    212,          /* lineNo */
    9,            /* colNo */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m" /* pName */
};

static emlrtBCInfo qg_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    218,          /* lineNo */
    9,            /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo rg_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    218,          /* lineNo */
    11,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo sg_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    219,          /* lineNo */
    9,            /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo tg_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    219,          /* lineNo */
    11,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtRTEInfo m_emlrtRTEI = {
    58,                   /* lineNo */
    23,                   /* colNo */
    "assertValidSizeArg", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2023a\\toolbox\\eml\\eml\\+coder\\+"
    "internal\\assertValidSizeArg.m" /* pName */
};

static emlrtBCInfo ug_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    177,          /* lineNo */
    34,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo vg_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    177,          /* lineNo */
    40,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo wg_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    177,          /* lineNo */
    54,           /* colNo */
    "muim",       /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo xg_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    177,          /* lineNo */
    11,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo yg_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    177,          /* lineNo */
    17,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo ah_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    178,          /* lineNo */
    34,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo bh_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    178,          /* lineNo */
    40,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo ch_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    178,          /* lineNo */
    54,           /* colNo */
    "muim",       /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo dh_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    178,          /* lineNo */
    11,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo eh_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    178,          /* lineNo */
    17,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo fh_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    160,          /* lineNo */
    16,           /* colNo */
    "xeta",       /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo gh_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    161,          /* lineNo */
    16,           /* colNo */
    "yeta",       /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo hh_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    186,          /* lineNo */
    34,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo ih_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    186,          /* lineNo */
    40,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo jh_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    186,          /* lineNo */
    56,           /* colNo */
    "muim",       /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo kh_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    186,          /* lineNo */
    11,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo lh_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    186,          /* lineNo */
    17,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo mh_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    187,          /* lineNo */
    34,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo nh_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    187,          /* lineNo */
    40,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo oh_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    187,          /* lineNo */
    56,           /* colNo */
    "muim",       /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo ph_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    187,          /* lineNo */
    11,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo qh_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    187,          /* lineNo */
    17,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo rh_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    199,          /* lineNo */
    32,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo sh_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    199,          /* lineNo */
    38,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo th_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    199,          /* lineNo */
    54,           /* colNo */
    "muim",       /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtDCInfo i_emlrtDCI = {
    199,          /* lineNo */
    54,           /* colNo */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    1                                                  /* checkKind */
};

static emlrtBCInfo uh_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    199,          /* lineNo */
    9,            /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo vh_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    199,          /* lineNo */
    15,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo wh_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    195,          /* lineNo */
    34,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo xh_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    195,          /* lineNo */
    40,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo yh_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    195,          /* lineNo */
    56,           /* colNo */
    "muim",       /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo ai_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    195,          /* lineNo */
    11,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo bi_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    195,          /* lineNo */
    17,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo ci_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    200,          /* lineNo */
    32,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo di_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    200,          /* lineNo */
    38,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo ei_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    200,          /* lineNo */
    54,           /* colNo */
    "muim",       /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo fi_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    200,          /* lineNo */
    9,            /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo gi_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    200,          /* lineNo */
    15,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo hi_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    196,          /* lineNo */
    34,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo ii_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    196,          /* lineNo */
    40,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo ji_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    196,          /* lineNo */
    56,           /* colNo */
    "muim",       /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo ki_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    196,          /* lineNo */
    11,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo li_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    196,          /* lineNo */
    17,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo mi_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    208,          /* lineNo */
    32,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo ni_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    208,          /* lineNo */
    38,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo oi_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    208,          /* lineNo */
    54,           /* colNo */
    "muim",       /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtDCInfo j_emlrtDCI = {
    208,          /* lineNo */
    54,           /* colNo */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    1                                                  /* checkKind */
};

static emlrtBCInfo pi_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    208,          /* lineNo */
    9,            /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo qi_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    208,          /* lineNo */
    15,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo ri_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    204,          /* lineNo */
    34,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo si_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    204,          /* lineNo */
    40,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo ti_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    204,          /* lineNo */
    56,           /* colNo */
    "muim",       /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo ui_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    204,          /* lineNo */
    11,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo vi_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    204,          /* lineNo */
    17,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo wi_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    209,          /* lineNo */
    32,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo xi_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    209,          /* lineNo */
    38,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo yi_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    209,          /* lineNo */
    54,           /* colNo */
    "muim",       /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo aj_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    209,          /* lineNo */
    9,            /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo bj_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    209,          /* lineNo */
    15,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo cj_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    205,          /* lineNo */
    34,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo dj_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    205,          /* lineNo */
    40,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo ej_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    205,          /* lineNo */
    56,           /* colNo */
    "muim",       /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo fj_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    205,          /* lineNo */
    11,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo gj_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    205,          /* lineNo */
    17,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo hj_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    213,          /* lineNo */
    34,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo ij_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    213,          /* lineNo */
    40,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo jj_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    213,          /* lineNo */
    54,           /* colNo */
    "muim",       /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo kj_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    213,          /* lineNo */
    11,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo lj_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    213,          /* lineNo */
    17,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo mj_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    214,          /* lineNo */
    34,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo nj_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    214,          /* lineNo */
    40,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo oj_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    214,          /* lineNo */
    54,           /* colNo */
    "muim",       /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo pj_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    214,          /* lineNo */
    11,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo qj_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    214,          /* lineNo */
    17,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo rj_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    220,          /* lineNo */
    9,            /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtDCInfo k_emlrtDCI = {
    220,          /* lineNo */
    9,            /* colNo */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    1                                                  /* checkKind */
};

static emlrtBCInfo sj_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    220,          /* lineNo */
    16,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtDCInfo l_emlrtDCI = {
    220,          /* lineNo */
    16,           /* colNo */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    1                                                  /* checkKind */
};

static emlrtBCInfo tj_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    221,          /* lineNo */
    9,            /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtDCInfo m_emlrtDCI = {
    221,          /* lineNo */
    9,            /* colNo */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    1                                                  /* checkKind */
};

static emlrtBCInfo uj_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    221,          /* lineNo */
    16,           /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtDCInfo n_emlrtDCI = {
    221,          /* lineNo */
    16,           /* colNo */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    1                                                  /* checkKind */
};

static emlrtBCInfo vj_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    167,          /* lineNo */
    3,            /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtBCInfo wj_emlrtBCI = {
    -1,           /* iFirst */
    -1,           /* iLast */
    166,          /* lineNo */
    3,            /* colNo */
    "lhs_tmp",    /* aName */
    "create_lhs", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m", /* pName */
    0                                                  /* checkKind */
};

static emlrtRTEInfo u_emlrtRTEI = {
    13,     /* lineNo */
    9,      /* colNo */
    "sqrt", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2023a\\toolbox\\eml\\lib\\matlab\\elfun\\sqrt.m" /* pName
                                                                       */
};

static emlrtRTEInfo v_emlrtRTEI = {
    20,                /* lineNo */
    1,                 /* colNo */
    "march_grid_fast", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m" /* pName */
};

static emlrtRTEInfo w_emlrtRTEI = {
    26,                /* lineNo */
    1,                 /* colNo */
    "march_grid_fast", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m" /* pName */
};

static emlrtRTEInfo x_emlrtRTEI = {
    33,                /* lineNo */
    1,                 /* colNo */
    "march_grid_fast", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m" /* pName */
};

static emlrtRTEInfo y_emlrtRTEI = {
    34,                /* lineNo */
    1,                 /* colNo */
    "march_grid_fast", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m" /* pName */
};

static emlrtRTEInfo ab_emlrtRTEI = {
    27,                /* lineNo */
    1,                 /* colNo */
    "march_grid_fast", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m" /* pName */
};

static emlrtRTEInfo bb_emlrtRTEI = {
    1,                 /* lineNo */
    35,                /* colNo */
    "march_grid_fast", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m" /* pName */
};

static emlrtRTEInfo cb_emlrtRTEI = {
    64,                /* lineNo */
    1,                 /* colNo */
    "march_grid_fast", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m" /* pName */
};

static emlrtRTEInfo db_emlrtRTEI = {
    65,                /* lineNo */
    1,                 /* colNo */
    "march_grid_fast", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m" /* pName */
};

static emlrtRTEInfo eb_emlrtRTEI = {
    66,                /* lineNo */
    1,                 /* colNo */
    "march_grid_fast", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m" /* pName */
};

static emlrtRTEInfo fb_emlrtRTEI = {
    67,                /* lineNo */
    1,                 /* colNo */
    "march_grid_fast", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m" /* pName */
};

static emlrtRTEInfo gb_emlrtRTEI = {
    151,               /* lineNo */
    1,                 /* colNo */
    "march_grid_fast", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m" /* pName */
};

static emlrtRTEInfo hb_emlrtRTEI = {
    223,               /* lineNo */
    1,                 /* colNo */
    "march_grid_fast", /* fName */
    "C:\\Users\\Will\\Documents\\MATLAB\\VT_Research\\2025\\grid_generation_"
    "stuff\\kt_airfoil_05_11_2026\\march_grid_fast.m" /* pName */
};

/* Function Declarations */
static void b_binary_expand_op(const emlrtStack *sp, emxArray_real_T *in1,
                               const emxArray_real_T *in2,
                               const emxArray_real_T *in3, const real_T in5[2]);

static void binary_expand_op(const emlrtStack *sp, emxArray_real_T *in1,
                             const emxArray_real_T *in2,
                             const emxArray_real_T *in3, const real_T in5[2]);

static int32_T
create_lhs(const emlrtStack *sp, real_T imax, const emxArray_real_T *xxi,
           const emxArray_real_T *yxi, const emxArray_real_T *xeta,
           const emxArray_real_T *yeta, const emxArray_real_T *alpha,
           const emxArray_real_T *muim, emxArray_real_T *lhs_d,
           emxArray_int32_T *lhs_colidx, emxArray_int32_T *lhs_rowidx,
           int32_T *lhs_n);

static void grid_metrics(const emlrtStack *sp, real_T imax,
                         const emxArray_real_T *x, const emxArray_real_T *y,
                         const emxArray_real_T *volume, emxArray_real_T *xxi,
                         emxArray_real_T *yxi, emxArray_real_T *xeta,
                         emxArray_real_T *yeta);

/* Function Definitions */
static void b_binary_expand_op(const emlrtStack *sp, emxArray_real_T *in1,
                               const emxArray_real_T *in2,
                               const emxArray_real_T *in3, const real_T in5[2])
{
  const real_T *in2_data;
  const real_T *in3_data;
  real_T *in1_data;
  int32_T i;
  int32_T in5_idx_0;
  int32_T loop_ub;
  int32_T stride_0_0;
  in3_data = in3->data;
  in2_data = in2->data;
  in5_idx_0 = (int32_T)in5[0];
  if (in5_idx_0 == 1) {
    loop_ub = in2->size[0];
  } else {
    loop_ub = in5_idx_0;
  }
  i = in1->size[0];
  in1->size[0] = loop_ub;
  emxEnsureCapacity_real_T(sp, in1, i, &x_emlrtRTEI);
  in1_data = in1->data;
  stride_0_0 = (in2->size[0] != 1);
  in5_idx_0 = (in5_idx_0 != 1);
  for (i = 0; i < loop_ub; i++) {
    in1_data[i] = in2_data[i * stride_0_0] + in3_data[2 * (i * in5_idx_0)];
  }
}

static void binary_expand_op(const emlrtStack *sp, emxArray_real_T *in1,
                             const emxArray_real_T *in2,
                             const emxArray_real_T *in3, const real_T in5[2])
{
  const real_T *in2_data;
  const real_T *in3_data;
  real_T *in1_data;
  int32_T i;
  int32_T in5_idx_0;
  int32_T loop_ub;
  int32_T stride_0_0;
  in3_data = in3->data;
  in2_data = in2->data;
  in5_idx_0 = (int32_T)in5[0];
  if (in5_idx_0 == 1) {
    loop_ub = in2->size[0];
  } else {
    loop_ub = in5_idx_0;
  }
  i = in1->size[0];
  in1->size[0] = loop_ub;
  emxEnsureCapacity_real_T(sp, in1, i, &y_emlrtRTEI);
  in1_data = in1->data;
  stride_0_0 = (in2->size[0] != 1);
  in5_idx_0 = (in5_idx_0 != 1);
  for (i = 0; i < loop_ub; i++) {
    in1_data[i] = in2_data[i * stride_0_0] + in3_data[2 * (i * in5_idx_0) + 1];
  }
}

static int32_T
create_lhs(const emlrtStack *sp, real_T imax, const emxArray_real_T *xxi,
           const emxArray_real_T *yxi, const emxArray_real_T *xeta,
           const emxArray_real_T *yeta, const emxArray_real_T *alpha,
           const emxArray_real_T *muim, emxArray_real_T *lhs_d,
           emxArray_int32_T *lhs_colidx, emxArray_int32_T *lhs_rowidx,
           int32_T *lhs_n)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack st;
  emxArray_real_T *lhs_tmp;
  const real_T *alpha_data;
  const real_T *muim_data;
  const real_T *xeta_data;
  const real_T *xxi_data;
  const real_T *yeta_data;
  const real_T *yxi_data;
  real_T BinvA_tmp;
  real_T fullk_tmp;
  real_T jacobian_tmp;
  real_T t;
  real_T xeta0;
  real_T *lhs_d_data;
  real_T *lhs_tmp_data;
  int32_T lhs_m;
  int32_T loop_ub;
  int32_T mInt;
  int32_T nInt;
  int32_T numalloc;
  int32_T row;
  int32_T *lhs_colidx_data;
  int32_T *lhs_rowidx_data;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  muim_data = muim->data;
  alpha_data = alpha->data;
  yeta_data = yeta->data;
  xeta_data = xeta->data;
  yxi_data = yxi->data;
  xxi_data = xxi->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  /* CREATE_LHS creates the LHS matrix for the standard central difference */
  /*  scheme plus fourth order implicit dissipation for smoothing */
  /*  Create square matrix with 1 on the diagonal */
  st.site = &q_emlrtRSI;
  fullk_tmp = 2.0 * imax;
  b_st.site = &t_emlrtRSI;
  if (fullk_tmp < 0.0) {
    t = 0.0;
  } else {
    t = fullk_tmp;
  }
  c_st.site = &v_emlrtRSI;
  if ((t != muDoubleScalarFloor(t)) || muDoubleScalarIsInf(t) ||
      (t > 2.147483647E+9)) {
    emlrtErrorWithMessageIdR2018a(
        &c_st, &m_emlrtRTEI, "Coder:MATLAB:NonIntegerInput",
        "Coder:MATLAB:NonIntegerInput", 4, 12, MIN_int32_T, 12, MAX_int32_T);
  }
  if (t <= 0.0) {
    jacobian_tmp = 0.0;
  } else {
    jacobian_tmp = t;
  }
  if (!(jacobian_tmp <= 2.147483647E+9)) {
    emlrtErrorWithMessageIdR2018a(&c_st, &n_emlrtRTEI, "Coder:MATLAB:pmaxsize",
                                  "Coder:MATLAB:pmaxsize", 0);
  }
  numalloc = (int32_T)t;
  emxInit_real_T(&st, &lhs_tmp, 2, &gb_emlrtRTEI);
  row = lhs_tmp->size[0] * lhs_tmp->size[1];
  lhs_tmp->size[0] = (int32_T)t;
  lhs_tmp->size[1] = (int32_T)t;
  emxEnsureCapacity_real_T(&st, lhs_tmp, row, &gb_emlrtRTEI);
  lhs_tmp_data = lhs_tmp->data;
  loop_ub = (int32_T)t * (int32_T)t;
  for (row = 0; row < loop_ub; row++) {
    lhs_tmp_data[row] = 0.0;
  }
  if ((int32_T)t > 0) {
    b_st.site = &u_emlrtRSI;
    if ((int32_T)t > 2147483646) {
      c_st.site = &w_emlrtRSI;
      check_forloop_overflow_error(&c_st);
    }
    for (loop_ub = 0; loop_ub < numalloc; loop_ub++) {
      lhs_tmp_data[loop_ub + lhs_tmp->size[0] * loop_ub] = 1.0;
    }
  }
  /*  Set LHS Based on Grid Type */
  /*  Apply Steger-Chaussee hyperbolic grid marching, eq 12 */
  /*  Kinsey and Barth pg 5 */
  row = (int32_T)((imax - 1.0) - 1.0);
  emlrtForLoopVectorCheckR2021a(2.0, 1.0, imax - 1.0, mxDOUBLE_CLASS,
                                (int32_T)((imax - 1.0) - 1.0), &g_emlrtRTEI,
                                (emlrtConstCTX)sp);
  for (loop_ub = 0; loop_ub < row; loop_ub++) {
    real_T BinvA_idx_0;
    real_T b_BinvA_tmp;
    real_T b_jacobian_tmp;
    boolean_T b;
    if ((loop_ub + 2 < 1) || (loop_ub + 2 > xxi->size[0])) {
      emlrtDynamicBoundsCheckR2012b(loop_ub + 2, 1, xxi->size[0], &sf_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    if ((loop_ub + 2 < 1) || (loop_ub + 2 > yxi->size[0])) {
      emlrtDynamicBoundsCheckR2012b(loop_ub + 2, 1, yxi->size[0], &tf_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    st.site = &r_emlrtRSI;
    b_st.site = &k_emlrtRSI;
    c_st.site = &l_emlrtRSI;
    st.site = &r_emlrtRSI;
    b_st.site = &k_emlrtRSI;
    c_st.site = &l_emlrtRSI;
    jacobian_tmp = xxi_data[loop_ub + 1];
    b_jacobian_tmp = yxi_data[loop_ub + 1];
    t = jacobian_tmp * jacobian_tmp + b_jacobian_tmp * b_jacobian_tmp;
    if (((int32_T)((uint32_T)loop_ub + 2U) < 1) ||
        ((int32_T)((uint32_T)loop_ub + 2U) > xeta->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)loop_ub + 2U), 1,
                                    xeta->size[0], &fh_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    xeta0 = xeta_data[loop_ub + 1] / t;
    if (((int32_T)((uint32_T)loop_ub + 2U) < 1) ||
        ((int32_T)((uint32_T)loop_ub + 2U) > yeta->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)loop_ub + 2U), 1,
                                    yeta->size[0], &gh_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    t = yeta_data[loop_ub + 1] / t;
    BinvA_tmp = b_jacobian_tmp * t;
    b_BinvA_tmp = jacobian_tmp * xeta0;
    BinvA_idx_0 = b_BinvA_tmp - BinvA_tmp;
    b_jacobian_tmp = jacobian_tmp * t + b_jacobian_tmp * xeta0;
    t = BinvA_tmp - b_BinvA_tmp;
    b = ((loop_ub + 2 < 1) || (loop_ub + 2 > alpha->size[0]));
    if (b) {
      emlrtDynamicBoundsCheckR2012b(loop_ub + 2, 1, alpha->size[0],
                                    &uf_emlrtBCI, (emlrtConstCTX)sp);
    }
    xeta0 = 2.0 * ((real_T)loop_ub + 2.0);
    if (((int32_T)(xeta0 - 1.0) < 1) ||
        ((int32_T)(xeta0 - 1.0) > lhs_tmp->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(xeta0 - 1.0), 1, lhs_tmp->size[0],
                                    &wj_emlrtBCI, (emlrtConstCTX)sp);
    }
    if (((int32_T)(xeta0 - 3.0) < 1) ||
        ((int32_T)(xeta0 - 3.0) > lhs_tmp->size[1])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(xeta0 - 3.0), 1, lhs_tmp->size[1],
                                    &wj_emlrtBCI, (emlrtConstCTX)sp);
    }
    jacobian_tmp = alpha_data[loop_ub + 1];
    lhs_tmp_data[((int32_T)(xeta0 - 1.0) +
                  lhs_tmp->size[0] * ((int32_T)(xeta0 - 3.0) - 1)) -
                 1] = -jacobian_tmp * BinvA_idx_0 / 2.0;
    if (((int32_T)((xeta0 - 1.0) + 1.0) < 1) ||
        ((int32_T)((xeta0 - 1.0) + 1.0) > lhs_tmp->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((xeta0 - 1.0) + 1.0), 1,
                                    lhs_tmp->size[0], &wj_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    if (((int32_T)(xeta0 - 3.0) < 1) ||
        ((int32_T)(xeta0 - 3.0) > lhs_tmp->size[1])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(xeta0 - 3.0), 1, lhs_tmp->size[1],
                                    &wj_emlrtBCI, (emlrtConstCTX)sp);
    }
    BinvA_tmp = -jacobian_tmp * b_jacobian_tmp / 2.0;
    lhs_tmp_data[((int32_T)((xeta0 - 1.0) + 1.0) +
                  lhs_tmp->size[0] * ((int32_T)(xeta0 - 3.0) - 1)) -
                 1] = BinvA_tmp;
    if (((int32_T)(xeta0 - 1.0) < 1) ||
        ((int32_T)(xeta0 - 1.0) > lhs_tmp->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(xeta0 - 1.0), 1, lhs_tmp->size[0],
                                    &wj_emlrtBCI, (emlrtConstCTX)sp);
    }
    if (((int32_T)((xeta0 - 3.0) + 1.0) < 1) ||
        ((int32_T)((xeta0 - 3.0) + 1.0) > lhs_tmp->size[1])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((xeta0 - 3.0) + 1.0), 1,
                                    lhs_tmp->size[1], &wj_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    lhs_tmp_data[((int32_T)(xeta0 - 1.0) +
                  lhs_tmp->size[0] * ((int32_T)((xeta0 - 3.0) + 1.0) - 1)) -
                 1] = BinvA_tmp;
    if (((int32_T)((xeta0 - 1.0) + 1.0) < 1) ||
        ((int32_T)((xeta0 - 1.0) + 1.0) > lhs_tmp->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((xeta0 - 1.0) + 1.0), 1,
                                    lhs_tmp->size[0], &wj_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    if (((int32_T)((xeta0 - 3.0) + 1.0) < 1) ||
        ((int32_T)((xeta0 - 3.0) + 1.0) > lhs_tmp->size[1])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((xeta0 - 3.0) + 1.0), 1,
                                    lhs_tmp->size[1], &wj_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    lhs_tmp_data[((int32_T)((xeta0 - 1.0) + 1.0) +
                  lhs_tmp->size[0] * ((int32_T)((xeta0 - 3.0) + 1.0) - 1)) -
                 1] = -jacobian_tmp * t / 2.0;
    if (((int32_T)(xeta0 - 1.0) < 1) ||
        ((int32_T)(xeta0 - 1.0) > lhs_tmp->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(xeta0 - 1.0), 1, lhs_tmp->size[0],
                                    &vj_emlrtBCI, (emlrtConstCTX)sp);
    }
    if (((int32_T)(xeta0 + 1.0) < 1) ||
        ((int32_T)(xeta0 + 1.0) > lhs_tmp->size[1])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(xeta0 + 1.0), 1, lhs_tmp->size[1],
                                    &vj_emlrtBCI, (emlrtConstCTX)sp);
    }
    lhs_tmp_data[((int32_T)(xeta0 - 1.0) +
                  lhs_tmp->size[0] * ((int32_T)(xeta0 + 1.0) - 1)) -
                 1] = jacobian_tmp * BinvA_idx_0 / 2.0;
    if (((int32_T)((xeta0 - 1.0) + 1.0) < 1) ||
        ((int32_T)((xeta0 - 1.0) + 1.0) > lhs_tmp->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((xeta0 - 1.0) + 1.0), 1,
                                    lhs_tmp->size[0], &vj_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    if (((int32_T)(xeta0 + 1.0) < 1) ||
        ((int32_T)(xeta0 + 1.0) > lhs_tmp->size[1])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(xeta0 + 1.0), 1, lhs_tmp->size[1],
                                    &vj_emlrtBCI, (emlrtConstCTX)sp);
    }
    BinvA_tmp = jacobian_tmp * b_jacobian_tmp / 2.0;
    lhs_tmp_data[((int32_T)((xeta0 - 1.0) + 1.0) +
                  lhs_tmp->size[0] * ((int32_T)(xeta0 + 1.0) - 1)) -
                 1] = BinvA_tmp;
    if (((int32_T)(xeta0 - 1.0) < 1) ||
        ((int32_T)(xeta0 - 1.0) > lhs_tmp->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(xeta0 - 1.0), 1, lhs_tmp->size[0],
                                    &vj_emlrtBCI, (emlrtConstCTX)sp);
    }
    if (((int32_T)((xeta0 + 1.0) + 1.0) < 1) ||
        ((int32_T)((xeta0 + 1.0) + 1.0) > lhs_tmp->size[1])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((xeta0 + 1.0) + 1.0), 1,
                                    lhs_tmp->size[1], &vj_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    lhs_tmp_data[((int32_T)(xeta0 - 1.0) +
                  lhs_tmp->size[0] * ((int32_T)((xeta0 + 1.0) + 1.0) - 1)) -
                 1] = BinvA_tmp;
    if (((int32_T)((xeta0 - 1.0) + 1.0) < 1) ||
        ((int32_T)((xeta0 - 1.0) + 1.0) > lhs_tmp->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((xeta0 - 1.0) + 1.0), 1,
                                    lhs_tmp->size[0], &vj_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    if (((int32_T)((xeta0 + 1.0) + 1.0) < 1) ||
        ((int32_T)((xeta0 + 1.0) + 1.0) > lhs_tmp->size[1])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((xeta0 + 1.0) + 1.0), 1,
                                    lhs_tmp->size[1], &vj_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    lhs_tmp_data[((int32_T)((xeta0 - 1.0) + 1.0) +
                  lhs_tmp->size[0] * ((int32_T)((xeta0 + 1.0) + 1.0) - 1)) -
                 1] = jacobian_tmp * t / 2.0;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b((emlrtConstCTX)sp);
    }
  }
  /*  Now add implicit smoothing to create pentadiagonal system */
  /*  Kinsey and Barth pg 5 */
  /*  Lowest block diagonal */
  row = (int32_T)((imax - 1.0) - 2.0);
  emlrtForLoopVectorCheckR2021a(3.0, 1.0, imax - 1.0, mxDOUBLE_CLASS,
                                (int32_T)((imax - 1.0) - 2.0), &h_emlrtRTEI,
                                (emlrtConstCTX)sp);
  for (loop_ub = 0; loop_ub < row; loop_ub++) {
    xeta0 = 2.0 * ((real_T)loop_ub + 3.0);
    if (((int32_T)(xeta0 - 1.0) < 1) ||
        ((int32_T)(xeta0 - 1.0) > lhs_tmp->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(xeta0 - 1.0), 1, lhs_tmp->size[0],
                                    &ug_emlrtBCI, (emlrtConstCTX)sp);
    }
    if (((int32_T)(xeta0 - 5.0) < 1) ||
        ((int32_T)(xeta0 - 5.0) > lhs_tmp->size[1])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(xeta0 - 5.0), 1, lhs_tmp->size[1],
                                    &vg_emlrtBCI, (emlrtConstCTX)sp);
    }
    if (((int32_T)((uint32_T)loop_ub + 3U) < 1) ||
        ((int32_T)((uint32_T)loop_ub + 3U) > muim->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)loop_ub + 3U), 1,
                                    muim->size[0], &wg_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    if (((int32_T)(xeta0 - 1.0) < 1) ||
        ((int32_T)(xeta0 - 1.0) > lhs_tmp->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(xeta0 - 1.0), 1, lhs_tmp->size[0],
                                    &xg_emlrtBCI, (emlrtConstCTX)sp);
    }
    if (((int32_T)(xeta0 - 5.0) < 1) ||
        ((int32_T)(xeta0 - 5.0) > lhs_tmp->size[1])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(xeta0 - 5.0), 1, lhs_tmp->size[1],
                                    &yg_emlrtBCI, (emlrtConstCTX)sp);
    }
    BinvA_tmp = muim_data[loop_ub + 2];
    lhs_tmp_data[((int32_T)(xeta0 - 1.0) +
                  lhs_tmp->size[0] * ((int32_T)(xeta0 - 5.0) - 1)) -
                 1] += BinvA_tmp;
    if (((int32_T)xeta0 < 1) || ((int32_T)xeta0 > lhs_tmp->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)xeta0, 1, lhs_tmp->size[0],
                                    &ah_emlrtBCI, (emlrtConstCTX)sp);
    }
    if (((int32_T)(xeta0 - 4.0) < 1) ||
        ((int32_T)(xeta0 - 4.0) > lhs_tmp->size[1])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(xeta0 - 4.0), 1, lhs_tmp->size[1],
                                    &bh_emlrtBCI, (emlrtConstCTX)sp);
    }
    if (((int32_T)((uint32_T)loop_ub + 3U) < 1) ||
        ((int32_T)((uint32_T)loop_ub + 3U) > muim->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)loop_ub + 3U), 1,
                                    muim->size[0], &ch_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    if (((int32_T)xeta0 < 1) || ((int32_T)xeta0 > lhs_tmp->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)xeta0, 1, lhs_tmp->size[0],
                                    &dh_emlrtBCI, (emlrtConstCTX)sp);
    }
    if (((int32_T)(xeta0 - 4.0) < 1) ||
        ((int32_T)(xeta0 - 4.0) > lhs_tmp->size[1])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(xeta0 - 4.0), 1, lhs_tmp->size[1],
                                    &eh_emlrtBCI, (emlrtConstCTX)sp);
    }
    lhs_tmp_data[((int32_T)xeta0 +
                  lhs_tmp->size[0] * ((int32_T)(xeta0 - 4.0) - 1)) -
                 1] += BinvA_tmp;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b((emlrtConstCTX)sp);
    }
  }
  /*  Lower block diagonal */
  if (lhs_tmp->size[0] < 3) {
    emlrtDynamicBoundsCheckR2012b(3, 1, lhs_tmp->size[0], &yf_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  if (lhs_tmp->size[1] < 1) {
    emlrtDynamicBoundsCheckR2012b(1, 1, lhs_tmp->size[1], &ag_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  if (lhs_tmp->size[0] < 3) {
    emlrtDynamicBoundsCheckR2012b(3, 1, lhs_tmp->size[0], &vf_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  if (lhs_tmp->size[1] < 1) {
    emlrtDynamicBoundsCheckR2012b(1, 1, lhs_tmp->size[1], &wf_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  if (muim->size[0] < 2) {
    emlrtDynamicBoundsCheckR2012b(2, 1, muim->size[0], &xf_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  lhs_tmp_data[2] -= 2.0 * muim_data[1];
  if (lhs_tmp->size[0] < 4) {
    emlrtDynamicBoundsCheckR2012b(4, 1, lhs_tmp->size[0], &eg_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  if (lhs_tmp->size[1] < 2) {
    emlrtDynamicBoundsCheckR2012b(2, 1, lhs_tmp->size[1], &fg_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  if (lhs_tmp->size[0] < 4) {
    emlrtDynamicBoundsCheckR2012b(4, 1, lhs_tmp->size[0], &bg_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  if (lhs_tmp->size[1] < 2) {
    emlrtDynamicBoundsCheckR2012b(2, 1, lhs_tmp->size[1], &cg_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  if (muim->size[0] < 2) {
    emlrtDynamicBoundsCheckR2012b(2, 1, muim->size[0], &dg_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  lhs_tmp_data[lhs_tmp->size[0] + 3] -= 2.0 * muim_data[1];
  emlrtForLoopVectorCheckR2021a(3.0, 1.0, imax - 1.0, mxDOUBLE_CLASS,
                                (int32_T)((imax - 1.0) - 2.0), &i_emlrtRTEI,
                                (emlrtConstCTX)sp);
  for (loop_ub = 0; loop_ub < row; loop_ub++) {
    xeta0 = 2.0 * ((real_T)loop_ub + 3.0);
    if (((int32_T)(xeta0 - 1.0) < 1) ||
        ((int32_T)(xeta0 - 1.0) > lhs_tmp->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(xeta0 - 1.0), 1, lhs_tmp->size[0],
                                    &hh_emlrtBCI, (emlrtConstCTX)sp);
    }
    if (((int32_T)(xeta0 - 3.0) < 1) ||
        ((int32_T)(xeta0 - 3.0) > lhs_tmp->size[1])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(xeta0 - 3.0), 1, lhs_tmp->size[1],
                                    &ih_emlrtBCI, (emlrtConstCTX)sp);
    }
    if (((int32_T)((uint32_T)loop_ub + 3U) < 1) ||
        ((int32_T)((uint32_T)loop_ub + 3U) > muim->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)loop_ub + 3U), 1,
                                    muim->size[0], &jh_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    if (((int32_T)(xeta0 - 1.0) < 1) ||
        ((int32_T)(xeta0 - 1.0) > lhs_tmp->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(xeta0 - 1.0), 1, lhs_tmp->size[0],
                                    &kh_emlrtBCI, (emlrtConstCTX)sp);
    }
    if (((int32_T)(xeta0 - 3.0) < 1) ||
        ((int32_T)(xeta0 - 3.0) > lhs_tmp->size[1])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(xeta0 - 3.0), 1, lhs_tmp->size[1],
                                    &lh_emlrtBCI, (emlrtConstCTX)sp);
    }
    BinvA_tmp = 4.0 * muim_data[loop_ub + 2];
    lhs_tmp_data[((int32_T)(xeta0 - 1.0) +
                  lhs_tmp->size[0] * ((int32_T)(xeta0 - 3.0) - 1)) -
                 1] -= BinvA_tmp;
    if (((int32_T)xeta0 < 1) || ((int32_T)xeta0 > lhs_tmp->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)xeta0, 1, lhs_tmp->size[0],
                                    &mh_emlrtBCI, (emlrtConstCTX)sp);
    }
    if (((int32_T)(xeta0 - 2.0) < 1) ||
        ((int32_T)(xeta0 - 2.0) > lhs_tmp->size[1])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(xeta0 - 2.0), 1, lhs_tmp->size[1],
                                    &nh_emlrtBCI, (emlrtConstCTX)sp);
    }
    if (((int32_T)((uint32_T)loop_ub + 3U) < 1) ||
        ((int32_T)((uint32_T)loop_ub + 3U) > muim->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)loop_ub + 3U), 1,
                                    muim->size[0], &oh_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    if (((int32_T)xeta0 < 1) || ((int32_T)xeta0 > lhs_tmp->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)xeta0, 1, lhs_tmp->size[0],
                                    &ph_emlrtBCI, (emlrtConstCTX)sp);
    }
    if (((int32_T)(xeta0 - 2.0) < 1) ||
        ((int32_T)(xeta0 - 2.0) > lhs_tmp->size[1])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(xeta0 - 2.0), 1, lhs_tmp->size[1],
                                    &qh_emlrtBCI, (emlrtConstCTX)sp);
    }
    lhs_tmp_data[((int32_T)xeta0 +
                  lhs_tmp->size[0] * ((int32_T)(xeta0 - 2.0) - 1)) -
                 1] -= BinvA_tmp;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b((emlrtConstCTX)sp);
    }
  }
  /*  Diagonal */
  if (lhs_tmp->size[0] < 3) {
    emlrtDynamicBoundsCheckR2012b(3, 1, lhs_tmp->size[0], &jg_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  if (lhs_tmp->size[1] < 3) {
    emlrtDynamicBoundsCheckR2012b(3, 1, lhs_tmp->size[1], &kg_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  if (lhs_tmp->size[0] < 3) {
    emlrtDynamicBoundsCheckR2012b(3, 1, lhs_tmp->size[0], &gg_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  if (lhs_tmp->size[1] < 3) {
    emlrtDynamicBoundsCheckR2012b(3, 1, lhs_tmp->size[1], &hg_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  if (muim->size[0] < 2) {
    emlrtDynamicBoundsCheckR2012b(2, 1, muim->size[0], &ig_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  lhs_tmp_data[lhs_tmp->size[0] * 2 + 2] += 5.0 * muim_data[1];
  if (lhs_tmp->size[0] < 4) {
    emlrtDynamicBoundsCheckR2012b(4, 1, lhs_tmp->size[0], &og_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  if (lhs_tmp->size[1] < 4) {
    emlrtDynamicBoundsCheckR2012b(4, 1, lhs_tmp->size[1], &pg_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  if (lhs_tmp->size[0] < 4) {
    emlrtDynamicBoundsCheckR2012b(4, 1, lhs_tmp->size[0], &lg_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  if (lhs_tmp->size[1] < 4) {
    emlrtDynamicBoundsCheckR2012b(4, 1, lhs_tmp->size[1], &mg_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  if (muim->size[0] < 2) {
    emlrtDynamicBoundsCheckR2012b(2, 1, muim->size[0], &ng_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  lhs_tmp_data[lhs_tmp->size[0] * 3 + 3] += 5.0 * muim_data[1];
  row = (int32_T)((imax - 2.0) - 2.0);
  emlrtForLoopVectorCheckR2021a(3.0, 1.0, imax - 2.0, mxDOUBLE_CLASS,
                                (int32_T)((imax - 2.0) - 2.0), &j_emlrtRTEI,
                                (emlrtConstCTX)sp);
  for (loop_ub = 0; loop_ub < row; loop_ub++) {
    xeta0 = 2.0 * ((real_T)loop_ub + 3.0);
    if (((int32_T)(xeta0 - 1.0) < 1) ||
        ((int32_T)(xeta0 - 1.0) > lhs_tmp->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(xeta0 - 1.0), 1, lhs_tmp->size[0],
                                    &wh_emlrtBCI, (emlrtConstCTX)sp);
    }
    if (((int32_T)(xeta0 - 1.0) < 1) ||
        ((int32_T)(xeta0 - 1.0) > lhs_tmp->size[1])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(xeta0 - 1.0), 1, lhs_tmp->size[1],
                                    &xh_emlrtBCI, (emlrtConstCTX)sp);
    }
    if (((int32_T)((uint32_T)loop_ub + 3U) < 1) ||
        ((int32_T)((uint32_T)loop_ub + 3U) > muim->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)loop_ub + 3U), 1,
                                    muim->size[0], &yh_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    if (((int32_T)(xeta0 - 1.0) < 1) ||
        ((int32_T)(xeta0 - 1.0) > lhs_tmp->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(xeta0 - 1.0), 1, lhs_tmp->size[0],
                                    &ai_emlrtBCI, (emlrtConstCTX)sp);
    }
    if (((int32_T)(xeta0 - 1.0) < 1) ||
        ((int32_T)(xeta0 - 1.0) > lhs_tmp->size[1])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(xeta0 - 1.0), 1, lhs_tmp->size[1],
                                    &bi_emlrtBCI, (emlrtConstCTX)sp);
    }
    BinvA_tmp = 6.0 * muim_data[loop_ub + 2];
    lhs_tmp_data[((int32_T)(xeta0 - 1.0) +
                  lhs_tmp->size[0] * ((int32_T)(xeta0 - 1.0) - 1)) -
                 1] += BinvA_tmp;
    if (((int32_T)xeta0 < 1) || ((int32_T)xeta0 > lhs_tmp->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)xeta0, 1, lhs_tmp->size[0],
                                    &hi_emlrtBCI, (emlrtConstCTX)sp);
    }
    if (((int32_T)xeta0 < 1) || ((int32_T)xeta0 > lhs_tmp->size[1])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)xeta0, 1, lhs_tmp->size[1],
                                    &ii_emlrtBCI, (emlrtConstCTX)sp);
    }
    if (((int32_T)((uint32_T)loop_ub + 3U) < 1) ||
        ((int32_T)((uint32_T)loop_ub + 3U) > muim->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)loop_ub + 3U), 1,
                                    muim->size[0], &ji_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    if (((int32_T)xeta0 < 1) || ((int32_T)xeta0 > lhs_tmp->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)xeta0, 1, lhs_tmp->size[0],
                                    &ki_emlrtBCI, (emlrtConstCTX)sp);
    }
    if (((int32_T)xeta0 < 1) || ((int32_T)xeta0 > lhs_tmp->size[1])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)xeta0, 1, lhs_tmp->size[1],
                                    &li_emlrtBCI, (emlrtConstCTX)sp);
    }
    lhs_tmp_data[((int32_T)xeta0 + lhs_tmp->size[0] * ((int32_T)xeta0 - 1)) -
                 1] += BinvA_tmp;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b((emlrtConstCTX)sp);
    }
  }
  xeta0 = 2.0 * (imax - 1.0);
  if (((int32_T)(xeta0 - 1.0) < 1) ||
      ((int32_T)(xeta0 - 1.0) > lhs_tmp->size[0])) {
    emlrtDynamicBoundsCheckR2012b((int32_T)(xeta0 - 1.0), 1, lhs_tmp->size[0],
                                  &rh_emlrtBCI, (emlrtConstCTX)sp);
  }
  if (((int32_T)(xeta0 - 1.0) < 1) ||
      ((int32_T)(xeta0 - 1.0) > lhs_tmp->size[1])) {
    emlrtDynamicBoundsCheckR2012b((int32_T)(xeta0 - 1.0), 1, lhs_tmp->size[1],
                                  &sh_emlrtBCI, (emlrtConstCTX)sp);
  }
  row = (int32_T)muDoubleScalarFloor(imax - 1.0);
  if (imax - 1.0 != row) {
    emlrtIntegerCheckR2012b(imax - 1.0, &i_emlrtDCI, (emlrtConstCTX)sp);
  }
  if (((int32_T)(imax - 1.0) < 1) || ((int32_T)(imax - 1.0) > muim->size[0])) {
    emlrtDynamicBoundsCheckR2012b((int32_T)(imax - 1.0), 1, muim->size[0],
                                  &th_emlrtBCI, (emlrtConstCTX)sp);
  }
  if (((int32_T)(xeta0 - 1.0) < 1) ||
      ((int32_T)(xeta0 - 1.0) > lhs_tmp->size[0])) {
    emlrtDynamicBoundsCheckR2012b((int32_T)(xeta0 - 1.0), 1, lhs_tmp->size[0],
                                  &uh_emlrtBCI, (emlrtConstCTX)sp);
  }
  if (((int32_T)(xeta0 - 1.0) < 1) ||
      ((int32_T)(xeta0 - 1.0) > lhs_tmp->size[1])) {
    emlrtDynamicBoundsCheckR2012b((int32_T)(xeta0 - 1.0), 1, lhs_tmp->size[1],
                                  &vh_emlrtBCI, (emlrtConstCTX)sp);
  }
  lhs_tmp_data[((int32_T)(xeta0 - 1.0) +
                lhs_tmp->size[0] * ((int32_T)(xeta0 - 1.0) - 1)) -
               1] += 5.0 * muim_data[(int32_T)(imax - 1.0) - 1];
  if (((int32_T)xeta0 < 1) || ((int32_T)xeta0 > lhs_tmp->size[0])) {
    emlrtDynamicBoundsCheckR2012b((int32_T)xeta0, 1, lhs_tmp->size[0],
                                  &ci_emlrtBCI, (emlrtConstCTX)sp);
  }
  if (((int32_T)xeta0 < 1) || ((int32_T)xeta0 > lhs_tmp->size[1])) {
    emlrtDynamicBoundsCheckR2012b((int32_T)xeta0, 1, lhs_tmp->size[1],
                                  &di_emlrtBCI, (emlrtConstCTX)sp);
  }
  if (((int32_T)(imax - 1.0) < 1) || ((int32_T)(imax - 1.0) > muim->size[0])) {
    emlrtDynamicBoundsCheckR2012b((int32_T)(imax - 1.0), 1, muim->size[0],
                                  &ei_emlrtBCI, (emlrtConstCTX)sp);
  }
  if (((int32_T)xeta0 < 1) || ((int32_T)xeta0 > lhs_tmp->size[0])) {
    emlrtDynamicBoundsCheckR2012b((int32_T)xeta0, 1, lhs_tmp->size[0],
                                  &fi_emlrtBCI, (emlrtConstCTX)sp);
  }
  if (((int32_T)xeta0 < 1) || ((int32_T)xeta0 > lhs_tmp->size[1])) {
    emlrtDynamicBoundsCheckR2012b((int32_T)xeta0, 1, lhs_tmp->size[1],
                                  &gi_emlrtBCI, (emlrtConstCTX)sp);
  }
  BinvA_tmp = muim_data[(int32_T)(imax - 1.0) - 1];
  lhs_tmp_data[((int32_T)xeta0 + lhs_tmp->size[0] * ((int32_T)xeta0 - 1)) -
               1] += 5.0 * BinvA_tmp;
  /*  Upper block diagonal */
  numalloc = (int32_T)((imax - 2.0) - 1.0);
  emlrtForLoopVectorCheckR2021a(2.0, 1.0, imax - 2.0, mxDOUBLE_CLASS,
                                (int32_T)((imax - 2.0) - 1.0), &k_emlrtRTEI,
                                (emlrtConstCTX)sp);
  for (loop_ub = 0; loop_ub < numalloc; loop_ub++) {
    t = 2.0 * ((real_T)loop_ub + 2.0);
    if (((int32_T)(t - 1.0) < 1) || ((int32_T)(t - 1.0) > lhs_tmp->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(t - 1.0), 1, lhs_tmp->size[0],
                                    &ri_emlrtBCI, (emlrtConstCTX)sp);
    }
    if (((int32_T)(t + 1.0) < 1) || ((int32_T)(t + 1.0) > lhs_tmp->size[1])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(t + 1.0), 1, lhs_tmp->size[1],
                                    &si_emlrtBCI, (emlrtConstCTX)sp);
    }
    if (((int32_T)((uint32_T)loop_ub + 2U) < 1) ||
        ((int32_T)((uint32_T)loop_ub + 2U) > muim->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)loop_ub + 2U), 1,
                                    muim->size[0], &ti_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    if (((int32_T)(t - 1.0) < 1) || ((int32_T)(t - 1.0) > lhs_tmp->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(t - 1.0), 1, lhs_tmp->size[0],
                                    &ui_emlrtBCI, (emlrtConstCTX)sp);
    }
    if (((int32_T)(t + 1.0) < 1) || ((int32_T)(t + 1.0) > lhs_tmp->size[1])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(t + 1.0), 1, lhs_tmp->size[1],
                                    &vi_emlrtBCI, (emlrtConstCTX)sp);
    }
    jacobian_tmp = 4.0 * muim_data[loop_ub + 1];
    lhs_tmp_data[((int32_T)(t - 1.0) +
                  lhs_tmp->size[0] * ((int32_T)(t + 1.0) - 1)) -
                 1] -= jacobian_tmp;
    if (((int32_T)t < 1) || ((int32_T)t > lhs_tmp->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)t, 1, lhs_tmp->size[0],
                                    &cj_emlrtBCI, (emlrtConstCTX)sp);
    }
    if (((int32_T)(t + 2.0) < 1) || ((int32_T)(t + 2.0) > lhs_tmp->size[1])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(t + 2.0), 1, lhs_tmp->size[1],
                                    &dj_emlrtBCI, (emlrtConstCTX)sp);
    }
    if (((int32_T)((uint32_T)loop_ub + 2U) < 1) ||
        ((int32_T)((uint32_T)loop_ub + 2U) > muim->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)loop_ub + 2U), 1,
                                    muim->size[0], &ej_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    if (((int32_T)t < 1) || ((int32_T)t > lhs_tmp->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)t, 1, lhs_tmp->size[0],
                                    &fj_emlrtBCI, (emlrtConstCTX)sp);
    }
    if (((int32_T)(t + 2.0) < 1) || ((int32_T)(t + 2.0) > lhs_tmp->size[1])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(t + 2.0), 1, lhs_tmp->size[1],
                                    &gj_emlrtBCI, (emlrtConstCTX)sp);
    }
    lhs_tmp_data[((int32_T)t + lhs_tmp->size[0] * ((int32_T)(t + 2.0) - 1)) -
                 1] -= jacobian_tmp;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b((emlrtConstCTX)sp);
    }
  }
  if (((int32_T)(xeta0 - 1.0) < 1) ||
      ((int32_T)(xeta0 - 1.0) > lhs_tmp->size[0])) {
    emlrtDynamicBoundsCheckR2012b((int32_T)(xeta0 - 1.0), 1, lhs_tmp->size[0],
                                  &mi_emlrtBCI, (emlrtConstCTX)sp);
  }
  if (((int32_T)(xeta0 + 1.0) < 1) ||
      ((int32_T)(xeta0 + 1.0) > lhs_tmp->size[1])) {
    emlrtDynamicBoundsCheckR2012b((int32_T)(xeta0 + 1.0), 1, lhs_tmp->size[1],
                                  &ni_emlrtBCI, (emlrtConstCTX)sp);
  }
  if (imax - 1.0 != row) {
    emlrtIntegerCheckR2012b(imax - 1.0, &j_emlrtDCI, (emlrtConstCTX)sp);
  }
  if (((int32_T)(imax - 1.0) < 1) || ((int32_T)(imax - 1.0) > muim->size[0])) {
    emlrtDynamicBoundsCheckR2012b((int32_T)(imax - 1.0), 1, muim->size[0],
                                  &oi_emlrtBCI, (emlrtConstCTX)sp);
  }
  if (((int32_T)(xeta0 - 1.0) < 1) ||
      ((int32_T)(xeta0 - 1.0) > lhs_tmp->size[0])) {
    emlrtDynamicBoundsCheckR2012b((int32_T)(xeta0 - 1.0), 1, lhs_tmp->size[0],
                                  &pi_emlrtBCI, (emlrtConstCTX)sp);
  }
  if (((int32_T)(xeta0 + 1.0) < 1) ||
      ((int32_T)(xeta0 + 1.0) > lhs_tmp->size[1])) {
    emlrtDynamicBoundsCheckR2012b((int32_T)(xeta0 + 1.0), 1, lhs_tmp->size[1],
                                  &qi_emlrtBCI, (emlrtConstCTX)sp);
  }
  lhs_tmp_data[((int32_T)(xeta0 - 1.0) +
                lhs_tmp->size[0] * ((int32_T)(xeta0 + 1.0) - 1)) -
               1] -= 2.0 * muim_data[(int32_T)(imax - 1.0) - 1];
  if (((int32_T)xeta0 < 1) || ((int32_T)xeta0 > lhs_tmp->size[0])) {
    emlrtDynamicBoundsCheckR2012b((int32_T)xeta0, 1, lhs_tmp->size[0],
                                  &wi_emlrtBCI, (emlrtConstCTX)sp);
  }
  if (((int32_T)(xeta0 + 2.0) < 1) ||
      ((int32_T)(xeta0 + 2.0) > lhs_tmp->size[1])) {
    emlrtDynamicBoundsCheckR2012b((int32_T)(xeta0 + 2.0), 1, lhs_tmp->size[1],
                                  &xi_emlrtBCI, (emlrtConstCTX)sp);
  }
  if (((int32_T)(imax - 1.0) < 1) || ((int32_T)(imax - 1.0) > muim->size[0])) {
    emlrtDynamicBoundsCheckR2012b((int32_T)(imax - 1.0), 1, muim->size[0],
                                  &yi_emlrtBCI, (emlrtConstCTX)sp);
  }
  if (((int32_T)xeta0 < 1) || ((int32_T)xeta0 > lhs_tmp->size[0])) {
    emlrtDynamicBoundsCheckR2012b((int32_T)xeta0, 1, lhs_tmp->size[0],
                                  &aj_emlrtBCI, (emlrtConstCTX)sp);
  }
  if (((int32_T)(xeta0 + 2.0) < 1) ||
      ((int32_T)(xeta0 + 2.0) > lhs_tmp->size[1])) {
    emlrtDynamicBoundsCheckR2012b((int32_T)(xeta0 + 2.0), 1, lhs_tmp->size[1],
                                  &bj_emlrtBCI, (emlrtConstCTX)sp);
  }
  lhs_tmp_data[((int32_T)xeta0 +
                lhs_tmp->size[0] * ((int32_T)(xeta0 + 2.0) - 1)) -
               1] -= 2.0 * BinvA_tmp;
  /*  Last diagonal */
  emlrtForLoopVectorCheckR2021a(2.0, 1.0, imax - 2.0, mxDOUBLE_CLASS,
                                (int32_T)((imax - 2.0) - 1.0), &l_emlrtRTEI,
                                (emlrtConstCTX)sp);
  for (loop_ub = 0; loop_ub < numalloc; loop_ub++) {
    xeta0 = 2.0 * ((real_T)loop_ub + 2.0);
    if (((int32_T)(xeta0 - 1.0) < 1) ||
        ((int32_T)(xeta0 - 1.0) > lhs_tmp->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(xeta0 - 1.0), 1, lhs_tmp->size[0],
                                    &hj_emlrtBCI, (emlrtConstCTX)sp);
    }
    if (((int32_T)(xeta0 + 3.0) < 1) ||
        ((int32_T)(xeta0 + 3.0) > lhs_tmp->size[1])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(xeta0 + 3.0), 1, lhs_tmp->size[1],
                                    &ij_emlrtBCI, (emlrtConstCTX)sp);
    }
    if (((int32_T)((uint32_T)loop_ub + 2U) < 1) ||
        ((int32_T)((uint32_T)loop_ub + 2U) > muim->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)loop_ub + 2U), 1,
                                    muim->size[0], &jj_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    if (((int32_T)(xeta0 - 1.0) < 1) ||
        ((int32_T)(xeta0 - 1.0) > lhs_tmp->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(xeta0 - 1.0), 1, lhs_tmp->size[0],
                                    &kj_emlrtBCI, (emlrtConstCTX)sp);
    }
    if (((int32_T)(xeta0 + 3.0) < 1) ||
        ((int32_T)(xeta0 + 3.0) > lhs_tmp->size[1])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(xeta0 + 3.0), 1, lhs_tmp->size[1],
                                    &lj_emlrtBCI, (emlrtConstCTX)sp);
    }
    BinvA_tmp = muim_data[loop_ub + 1];
    lhs_tmp_data[((int32_T)(xeta0 - 1.0) +
                  lhs_tmp->size[0] * ((int32_T)(xeta0 + 3.0) - 1)) -
                 1] += BinvA_tmp;
    if (((int32_T)xeta0 < 1) || ((int32_T)xeta0 > lhs_tmp->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)xeta0, 1, lhs_tmp->size[0],
                                    &mj_emlrtBCI, (emlrtConstCTX)sp);
    }
    if (((int32_T)(xeta0 + 4.0) < 1) ||
        ((int32_T)(xeta0 + 4.0) > lhs_tmp->size[1])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(xeta0 + 4.0), 1, lhs_tmp->size[1],
                                    &nj_emlrtBCI, (emlrtConstCTX)sp);
    }
    if (((int32_T)((uint32_T)loop_ub + 2U) < 1) ||
        ((int32_T)((uint32_T)loop_ub + 2U) > muim->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)loop_ub + 2U), 1,
                                    muim->size[0], &oj_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    if (((int32_T)xeta0 < 1) || ((int32_T)xeta0 > lhs_tmp->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)xeta0, 1, lhs_tmp->size[0],
                                    &pj_emlrtBCI, (emlrtConstCTX)sp);
    }
    if (((int32_T)(xeta0 + 4.0) < 1) ||
        ((int32_T)(xeta0 + 4.0) > lhs_tmp->size[1])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(xeta0 + 4.0), 1, lhs_tmp->size[1],
                                    &qj_emlrtBCI, (emlrtConstCTX)sp);
    }
    lhs_tmp_data[((int32_T)xeta0 +
                  lhs_tmp->size[0] * ((int32_T)(xeta0 + 4.0) - 1)) -
                 1] += BinvA_tmp;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b((emlrtConstCTX)sp);
    }
  }
  /*  Boundary Condition */
  if (lhs_tmp->size[0] < 2) {
    emlrtDynamicBoundsCheckR2012b(2, 1, lhs_tmp->size[0], &qg_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  if (lhs_tmp->size[1] < 4) {
    emlrtDynamicBoundsCheckR2012b(4, 1, lhs_tmp->size[1], &rg_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  lhs_tmp_data[lhs_tmp->size[0] * 3 + 1] = -2.0;
  if (lhs_tmp->size[0] < 2) {
    emlrtDynamicBoundsCheckR2012b(2, 1, lhs_tmp->size[0], &sg_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  if (lhs_tmp->size[1] < 6) {
    emlrtDynamicBoundsCheckR2012b(6, 1, lhs_tmp->size[1], &tg_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  lhs_tmp_data[lhs_tmp->size[0] * 5 + 1] = 1.0;
  row = (int32_T)muDoubleScalarFloor(fullk_tmp);
  if (fullk_tmp != row) {
    emlrtIntegerCheckR2012b(fullk_tmp, &k_emlrtDCI, (emlrtConstCTX)sp);
  }
  numalloc = (int32_T)fullk_tmp;
  if ((fullk_tmp < 1.0) || (numalloc > lhs_tmp->size[0])) {
    emlrtDynamicBoundsCheckR2012b((int32_T)fullk_tmp, 1, lhs_tmp->size[0],
                                  &rj_emlrtBCI, (emlrtConstCTX)sp);
  }
  xeta0 = 2.0 * imax - 2.0;
  if (xeta0 != (int32_T)muDoubleScalarFloor(xeta0)) {
    emlrtIntegerCheckR2012b(xeta0, &l_emlrtDCI, (emlrtConstCTX)sp);
  }
  if (((int32_T)xeta0 < 1) || ((int32_T)xeta0 > lhs_tmp->size[1])) {
    emlrtDynamicBoundsCheckR2012b((int32_T)xeta0, 1, lhs_tmp->size[1],
                                  &sj_emlrtBCI, (emlrtConstCTX)sp);
  }
  lhs_tmp_data[(numalloc + lhs_tmp->size[0] * ((int32_T)xeta0 - 1)) - 1] = -2.0;
  if (numalloc != row) {
    emlrtIntegerCheckR2012b(fullk_tmp, &m_emlrtDCI, (emlrtConstCTX)sp);
  }
  if ((fullk_tmp < 1.0) || (numalloc > lhs_tmp->size[0])) {
    emlrtDynamicBoundsCheckR2012b((int32_T)fullk_tmp, 1, lhs_tmp->size[0],
                                  &tj_emlrtBCI, (emlrtConstCTX)sp);
  }
  xeta0 = 2.0 * imax - 4.0;
  if (xeta0 != (int32_T)muDoubleScalarFloor(xeta0)) {
    emlrtIntegerCheckR2012b(xeta0, &n_emlrtDCI, (emlrtConstCTX)sp);
  }
  if (((int32_T)xeta0 < 1) || ((int32_T)xeta0 > lhs_tmp->size[1])) {
    emlrtDynamicBoundsCheckR2012b((int32_T)xeta0, 1, lhs_tmp->size[1],
                                  &uj_emlrtBCI, (emlrtConstCTX)sp);
  }
  lhs_tmp_data[(numalloc + lhs_tmp->size[0] * ((int32_T)xeta0 - 1)) - 1] = 1.0;
  st.site = &s_emlrtRSI;
  b_st.site = &x_emlrtRSI;
  mInt = lhs_tmp->size[0];
  nInt = lhs_tmp->size[1];
  c_st.site = &y_emlrtRSI;
  if (lhs_tmp->size[0] >= MAX_int32_T) {
    emlrtErrorWithMessageIdR2018a(
        &c_st, &f_emlrtRTEI, "Coder:toolbox:SparseMaxSize",
        "Coder:toolbox:SparseMaxSize", 2, 12, MAX_int32_T);
  }
  c_st.site = &ab_emlrtRSI;
  if (lhs_tmp->size[1] >= MAX_int32_T) {
    emlrtErrorWithMessageIdR2018a(
        &c_st, &f_emlrtRTEI, "Coder:toolbox:SparseMaxSize",
        "Coder:toolbox:SparseMaxSize", 2, 12, MAX_int32_T);
  }
  numalloc = 0;
  row = lhs_tmp->size[0] * lhs_tmp->size[1];
  for (loop_ub = 0; loop_ub < row; loop_ub++) {
    if (lhs_tmp_data[loop_ub] != 0.0) {
      numalloc++;
    }
  }
  lhs_m = lhs_tmp->size[0];
  *lhs_n = lhs_tmp->size[1];
  numalloc = muIntScalarMax_sint32(numalloc, 1);
  row = lhs_d->size[0];
  lhs_d->size[0] = numalloc;
  emxEnsureCapacity_real_T(&b_st, lhs_d, row, &hb_emlrtRTEI);
  lhs_d_data = lhs_d->data;
  for (row = 0; row < numalloc; row++) {
    lhs_d_data[row] = 0.0;
  }
  row = lhs_colidx->size[0];
  lhs_colidx->size[0] = lhs_tmp->size[1] + 1;
  emxEnsureCapacity_int32_T(&b_st, lhs_colidx, row, &hb_emlrtRTEI);
  lhs_colidx_data = lhs_colidx->data;
  loop_ub = lhs_tmp->size[1];
  for (row = 0; row <= loop_ub; row++) {
    lhs_colidx_data[row] = 0;
  }
  lhs_colidx_data[0] = 1;
  row = lhs_rowidx->size[0];
  lhs_rowidx->size[0] = numalloc;
  emxEnsureCapacity_int32_T(&b_st, lhs_rowidx, row, &hb_emlrtRTEI);
  lhs_rowidx_data = lhs_rowidx->data;
  for (row = 0; row < numalloc; row++) {
    lhs_rowidx_data[row] = 0;
  }
  lhs_rowidx_data[0] = 1;
  numalloc = 0;
  c_st.site = &bb_emlrtRSI;
  for (loop_ub = 0; loop_ub < nInt; loop_ub++) {
    c_st.site = &cb_emlrtRSI;
    if (mInt > 2147483646) {
      d_st.site = &w_emlrtRSI;
      check_forloop_overflow_error(&d_st);
    }
    for (row = 0; row < mInt; row++) {
      t = lhs_tmp_data[row + lhs_tmp->size[0] * loop_ub];
      if (t != 0.0) {
        lhs_rowidx_data[numalloc] = row + 1;
        lhs_d_data[numalloc] = t;
        numalloc++;
      }
    }
    lhs_colidx_data[loop_ub + 1] = numalloc + 1;
  }
  emxFree_real_T(&b_st, &lhs_tmp);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
  return lhs_m;
}

static void grid_metrics(const emlrtStack *sp, real_T imax,
                         const emxArray_real_T *x, const emxArray_real_T *y,
                         const emxArray_real_T *volume, emxArray_real_T *xxi,
                         emxArray_real_T *yxi, emxArray_real_T *xeta,
                         emxArray_real_T *yeta)
{
  emlrtStack st;
  const real_T *volume_data;
  const real_T *x_data;
  const real_T *y_data;
  real_T b_factor_tmp;
  real_T factor;
  real_T factor_tmp;
  real_T xeta_tmp;
  real_T *xeta_data;
  real_T *xxi_data;
  real_T *yeta_data;
  real_T *yxi_data;
  int32_T b_i;
  int32_T i;
  int32_T loop_ub_tmp_tmp;
  boolean_T b;
  st.prev = sp;
  st.tls = sp->tls;
  volume_data = volume->data;
  y_data = y->data;
  x_data = x->data;
  /* GRID_METRICS compute the current xi metrics and then estimate the eta */
  /*  metrics for the hyperbolic marching scheme */
  if (!(imax >= 0.0)) {
    emlrtNonNegativeCheckR2012b(imax, &e_emlrtDCI, (emlrtConstCTX)sp);
  }
  i = (int32_T)muDoubleScalarFloor(imax);
  if (imax != i) {
    emlrtIntegerCheckR2012b(imax, &d_emlrtDCI, (emlrtConstCTX)sp);
  }
  b_i = xxi->size[0];
  xxi->size[0] = (int32_T)imax;
  emxEnsureCapacity_real_T(sp, xxi, b_i, &cb_emlrtRTEI);
  xxi_data = xxi->data;
  if (imax != i) {
    emlrtIntegerCheckR2012b(imax, &d_emlrtDCI, (emlrtConstCTX)sp);
  }
  loop_ub_tmp_tmp = (int32_T)imax;
  for (b_i = 0; b_i < loop_ub_tmp_tmp; b_i++) {
    xxi_data[b_i] = 0.0;
  }
  if (loop_ub_tmp_tmp != i) {
    emlrtIntegerCheckR2012b(imax, &f_emlrtDCI, (emlrtConstCTX)sp);
  }
  b_i = yxi->size[0];
  yxi->size[0] = loop_ub_tmp_tmp;
  emxEnsureCapacity_real_T(sp, yxi, b_i, &db_emlrtRTEI);
  yxi_data = yxi->data;
  if (loop_ub_tmp_tmp != i) {
    emlrtIntegerCheckR2012b(imax, &f_emlrtDCI, (emlrtConstCTX)sp);
  }
  for (b_i = 0; b_i < loop_ub_tmp_tmp; b_i++) {
    yxi_data[b_i] = 0.0;
  }
  if (loop_ub_tmp_tmp != i) {
    emlrtIntegerCheckR2012b(imax, &g_emlrtDCI, (emlrtConstCTX)sp);
  }
  b_i = xeta->size[0];
  xeta->size[0] = loop_ub_tmp_tmp;
  emxEnsureCapacity_real_T(sp, xeta, b_i, &eb_emlrtRTEI);
  xeta_data = xeta->data;
  if (loop_ub_tmp_tmp != i) {
    emlrtIntegerCheckR2012b(imax, &g_emlrtDCI, (emlrtConstCTX)sp);
  }
  for (b_i = 0; b_i < loop_ub_tmp_tmp; b_i++) {
    xeta_data[b_i] = 0.0;
  }
  if (loop_ub_tmp_tmp != i) {
    emlrtIntegerCheckR2012b(imax, &h_emlrtDCI, (emlrtConstCTX)sp);
  }
  b_i = yeta->size[0];
  yeta->size[0] = loop_ub_tmp_tmp;
  emxEnsureCapacity_real_T(sp, yeta, b_i, &fb_emlrtRTEI);
  yeta_data = yeta->data;
  if (loop_ub_tmp_tmp != i) {
    emlrtIntegerCheckR2012b(imax, &h_emlrtDCI, (emlrtConstCTX)sp);
  }
  for (i = 0; i < loop_ub_tmp_tmp; i++) {
    yeta_data[i] = 0.0;
  }
  i = loop_ub_tmp_tmp - 2;
  emlrtForLoopVectorCheckR2021a(2.0, 1.0, imax - 1.0, mxDOUBLE_CLASS,
                                (int32_T)imax - 2, &e_emlrtRTEI,
                                (emlrtConstCTX)sp);
  for (b_i = 0; b_i < i; b_i++) {
    /*  xi metrics are set by the previous grid line */
    if ((int32_T)((uint32_T)b_i + 3U) > x->size[0]) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)b_i + 3U), 1,
                                    x->size[0], &te_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    if ((int32_T)((uint32_T)b_i + 1U) > x->size[0]) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)b_i + 1U), 1,
                                    x->size[0], &ue_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    if ((int32_T)((uint32_T)b_i + 2U) > xxi->size[0]) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)b_i + 2U), 1,
                                    xxi->size[0], &ve_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    xxi_data[b_i + 1] = 0.5 * (x_data[b_i + 2] - x_data[b_i]);
    if ((int32_T)((uint32_T)b_i + 3U) > y->size[0]) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)b_i + 3U), 1,
                                    y->size[0], &we_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    if ((int32_T)((uint32_T)b_i + 1U) > y->size[0]) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)b_i + 1U), 1,
                                    y->size[0], &xe_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    if ((int32_T)((uint32_T)b_i + 2U) > yxi->size[0]) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)b_i + 2U), 1,
                                    yxi->size[0], &ye_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    yxi_data[b_i + 1] = 0.5 * (y_data[b_i + 2] - y_data[b_i]);
    /*  eta metrics are guessed from the xi metrics and volumes */
    st.site = &m_emlrtRSI;
    if (b_i + 2 > xxi->size[0]) {
      emlrtDynamicBoundsCheckR2012b(b_i + 2, 1, xxi->size[0], &yd_emlrtBCI,
                                    &st);
    }
    st.site = &m_emlrtRSI;
    if (b_i + 2 > yxi->size[0]) {
      emlrtDynamicBoundsCheckR2012b(b_i + 2, 1, yxi->size[0], &ae_emlrtBCI,
                                    &st);
    }
    factor_tmp = xxi_data[b_i + 1];
    b_factor_tmp = yxi_data[b_i + 1];
    factor = factor_tmp * factor_tmp + b_factor_tmp * b_factor_tmp;
    if ((int32_T)((uint32_T)b_i + 2U) > yxi->size[0]) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)b_i + 2U), 1,
                                    yxi->size[0], &af_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    if ((int32_T)((uint32_T)b_i + 2U) > volume->size[0]) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)b_i + 2U), 1,
                                    volume->size[0], &bf_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    if ((int32_T)((uint32_T)b_i + 2U) > xeta->size[0]) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)b_i + 2U), 1,
                                    xeta->size[0], &cf_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    xeta_tmp = volume_data[b_i + 1];
    xeta_data[b_i + 1] = -yxi_data[b_i + 1] * xeta_tmp / factor;
    if ((int32_T)((uint32_T)b_i + 2U) > xxi->size[0]) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)b_i + 2U), 1,
                                    xxi->size[0], &df_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    if ((int32_T)((uint32_T)b_i + 2U) > volume->size[0]) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)b_i + 2U), 1,
                                    volume->size[0], &ef_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    if ((int32_T)((uint32_T)b_i + 2U) > yeta->size[0]) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)b_i + 2U), 1,
                                    yeta->size[0], &ff_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    yeta_data[b_i + 1] = xxi_data[b_i + 1] * xeta_tmp / factor;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b((emlrtConstCTX)sp);
    }
  }
  if (xxi->size[0] < 1) {
    emlrtDynamicBoundsCheckR2012b(1, 1, xxi->size[0], &ee_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  if (x->size[0] < 1) {
    emlrtDynamicBoundsCheckR2012b(1, 1, x->size[0], &be_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  if (x->size[0] < 2) {
    emlrtDynamicBoundsCheckR2012b(2, 1, x->size[0], &ce_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  if (x->size[0] < 3) {
    emlrtDynamicBoundsCheckR2012b(3, 1, x->size[0], &de_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  xxi_data[0] = 0.5 * ((-3.0 * x_data[0] + 4.0 * x_data[1]) - x_data[2]);
  if (yxi->size[0] < 1) {
    emlrtDynamicBoundsCheckR2012b(1, 1, yxi->size[0], &ie_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  if (y->size[0] < 1) {
    emlrtDynamicBoundsCheckR2012b(1, 1, y->size[0], &fe_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  if (y->size[0] < 2) {
    emlrtDynamicBoundsCheckR2012b(2, 1, y->size[0], &ge_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  if (y->size[0] < 3) {
    emlrtDynamicBoundsCheckR2012b(3, 1, y->size[0], &he_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  yxi_data[0] = 0.5 * ((-3.0 * y_data[0] + 4.0 * y_data[1]) - y_data[2]);
  if (xxi->size[0] < 1) {
    emlrtDynamicBoundsCheckR2012b(1, 1, xxi->size[0], &je_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  if (yxi->size[0] < 1) {
    emlrtDynamicBoundsCheckR2012b(1, 1, yxi->size[0], &ke_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  factor = xxi_data[0] * xxi_data[0] + yxi_data[0] * yxi_data[0];
  if (xeta->size[0] < 1) {
    emlrtDynamicBoundsCheckR2012b(1, 1, xeta->size[0], &re_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  if (yxi->size[0] < 1) {
    emlrtDynamicBoundsCheckR2012b(1, 1, yxi->size[0], &le_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  if (volume->size[0] < 1) {
    emlrtDynamicBoundsCheckR2012b(1, 1, volume->size[0], &me_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  xeta_data[0] = -yxi_data[0] * volume_data[0] / factor;
  if (yeta->size[0] < 1) {
    emlrtDynamicBoundsCheckR2012b(1, 1, yeta->size[0], &se_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  if (xxi->size[0] < 1) {
    emlrtDynamicBoundsCheckR2012b(1, 1, xxi->size[0], &ne_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  if (volume->size[0] < 1) {
    emlrtDynamicBoundsCheckR2012b(1, 1, volume->size[0], &oe_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  yeta_data[0] = xxi_data[0] * volume_data[0] / factor;
  if ((imax < 1.0) || (loop_ub_tmp_tmp > x->size[0])) {
    emlrtDynamicBoundsCheckR2012b((int32_T)imax, 1, x->size[0], &gf_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  if ((imax - 1.0 < 1.0) || (loop_ub_tmp_tmp - 1 > x->size[0])) {
    emlrtDynamicBoundsCheckR2012b((int32_T)imax - 1, 1, x->size[0],
                                  &hf_emlrtBCI, (emlrtConstCTX)sp);
  }
  if ((imax - 2.0 < 1.0) || (loop_ub_tmp_tmp - 2 > x->size[0])) {
    emlrtDynamicBoundsCheckR2012b((int32_T)imax - 2, 1, x->size[0],
                                  &if_emlrtBCI, (emlrtConstCTX)sp);
  }
  if ((imax < 1.0) || (loop_ub_tmp_tmp > xxi->size[0])) {
    emlrtDynamicBoundsCheckR2012b((int32_T)imax, 1, xxi->size[0], &jf_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  xxi_data[loop_ub_tmp_tmp - 1] = -0.5 * ((-3.0 * x_data[loop_ub_tmp_tmp - 1] +
                                           4.0 * x_data[loop_ub_tmp_tmp - 2]) -
                                          x_data[loop_ub_tmp_tmp - 3]);
  if ((imax < 1.0) || (loop_ub_tmp_tmp > y->size[0])) {
    emlrtDynamicBoundsCheckR2012b((int32_T)imax, 1, y->size[0], &kf_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  if ((imax - 1.0 < 1.0) || (loop_ub_tmp_tmp - 1 > y->size[0])) {
    emlrtDynamicBoundsCheckR2012b((int32_T)imax - 1, 1, y->size[0],
                                  &lf_emlrtBCI, (emlrtConstCTX)sp);
  }
  if ((imax - 2.0 < 1.0) || (loop_ub_tmp_tmp - 2 > y->size[0])) {
    emlrtDynamicBoundsCheckR2012b((int32_T)imax - 2, 1, y->size[0],
                                  &mf_emlrtBCI, (emlrtConstCTX)sp);
  }
  if ((imax < 1.0) || (loop_ub_tmp_tmp > yxi->size[0])) {
    emlrtDynamicBoundsCheckR2012b((int32_T)imax, 1, yxi->size[0], &nf_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  yxi_data[loop_ub_tmp_tmp - 1] = -0.5 * ((-3.0 * y_data[loop_ub_tmp_tmp - 1] +
                                           4.0 * y_data[loop_ub_tmp_tmp - 2]) -
                                          y_data[loop_ub_tmp_tmp - 3]);
  st.site = &o_emlrtRSI;
  b = ((imax < 1.0) || (loop_ub_tmp_tmp > xxi->size[0]));
  if (b) {
    emlrtDynamicBoundsCheckR2012b((int32_T)imax, 1, xxi->size[0], &pe_emlrtBCI,
                                  &st);
  }
  st.site = &o_emlrtRSI;
  b = ((imax < 1.0) || (loop_ub_tmp_tmp > yxi->size[0]));
  if (b) {
    emlrtDynamicBoundsCheckR2012b((int32_T)imax, 1, yxi->size[0], &qe_emlrtBCI,
                                  &st);
  }
  factor_tmp = xxi_data[loop_ub_tmp_tmp - 1];
  b_factor_tmp = yxi_data[loop_ub_tmp_tmp - 1];
  factor = factor_tmp * factor_tmp + b_factor_tmp * b_factor_tmp;
  if ((imax < 1.0) || (loop_ub_tmp_tmp > volume->size[0])) {
    emlrtDynamicBoundsCheckR2012b((int32_T)imax, 1, volume->size[0],
                                  &of_emlrtBCI, (emlrtConstCTX)sp);
  }
  if ((imax < 1.0) || (loop_ub_tmp_tmp > xeta->size[0])) {
    emlrtDynamicBoundsCheckR2012b((int32_T)imax, 1, xeta->size[0], &pf_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  xeta_tmp = volume_data[loop_ub_tmp_tmp - 1];
  xeta_data[loop_ub_tmp_tmp - 1] = -b_factor_tmp * xeta_tmp / factor;
  if ((imax < 1.0) || (loop_ub_tmp_tmp > volume->size[0])) {
    emlrtDynamicBoundsCheckR2012b((int32_T)imax, 1, volume->size[0],
                                  &qf_emlrtBCI, (emlrtConstCTX)sp);
  }
  if ((imax < 1.0) || (loop_ub_tmp_tmp > yeta->size[0])) {
    emlrtDynamicBoundsCheckR2012b((int32_T)imax, 1, yeta->size[0], &rf_emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  yeta_data[loop_ub_tmp_tmp - 1] = factor_tmp * xeta_tmp / factor;
}

void march_grid_fast(const emlrtStack *sp, real_T imax, real_T jmax,
                     const emxArray_real_T *x, const emxArray_real_T *y,
                     const emxArray_real_T *mu, const emxArray_real_T *muim,
                     const emxArray_real_T *alpham,
                     const emxArray_real_T *scale, const emxArray_real_T *rj,
                     const emxArray_real_T *rjm1, emxArray_real_T *x_update,
                     emxArray_real_T *y_update)
{
  emlrtStack b_st;
  emlrtStack st;
  emxArray_int32_T *lhs_colidx;
  emxArray_int32_T *lhs_rowidx;
  emxArray_real_T *lhs_d;
  emxArray_real_T *rhs;
  emxArray_real_T *volume;
  emxArray_real_T *xeta;
  emxArray_real_T *xxi;
  emxArray_real_T *yeta;
  emxArray_real_T *yxi;
  real_T dv[2];
  const real_T *mu_data;
  const real_T *rj_data;
  const real_T *rjm1_data;
  const real_T *scale_data;
  const real_T *x_data;
  const real_T *y_data;
  real_T arc_length;
  real_T dl;
  real_T vblend;
  real_T *rhs_data;
  real_T *volume_data;
  real_T *x_update_data;
  real_T *xeta_data;
  real_T *xxi_data;
  real_T *yxi_data;
  int32_T i;
  int32_T loop_ub_tmp;
  int32_T n;
  int32_T nx;
  uint32_T u;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  rjm1_data = rjm1->data;
  rj_data = rj->data;
  scale_data = scale->data;
  mu_data = mu->data;
  y_data = y->data;
  x_data = x->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  /* MARCH_GRID Main driver of grid marching in the eta direction */
  /*  INPUTS */
  /*  j       : current j index             (integer) */
  /*  jmax    : maximum j index             (integer) */
  /*  x       : x coordinates at level j-1  (double size [1,imax]) */
  /*  y       : y coordinates at level j-1  (double size [1,imax]) */
  /*  mu      : explicit scheme dissipation (double size [1,imax]) */
  /*  muim    : implicit scheme dissipation (double size [1,imax]) */
  /*  alpham  : implicit scheme parameter   (double size [1,imax]) */
  /*           (alpham = 0.5 ~ trapezoid rule) */
  /*           (alpham = 1.0 ~ backward Euler) */
  /*  scale   : volume scaling factor       (double size [1,imax]) */
  /*  ds      : radial spacing              (double size [1,imax]) */
  /*  Now start building the next grid line */
  /*  Calculate the estimated grid volumes */
  st.site = &emlrtRSI;
  /* VOLUME Based on scaling and the local arc length, calculate the estimated
   */
  /*  volumes, V^{0} */
  emxInit_real_T(&st, &volume, 1, &v_emlrtRTEI);
  if (!(imax >= 0.0)) {
    emlrtNonNegativeCheckR2012b(imax, &b_emlrtDCI, &st);
  }
  n = (int32_T)muDoubleScalarFloor(imax);
  if (imax != n) {
    emlrtIntegerCheckR2012b(imax, &c_emlrtDCI, &st);
  }
  loop_ub_tmp = (int32_T)imax;
  nx = volume->size[0];
  volume->size[0] = loop_ub_tmp;
  emxEnsureCapacity_real_T(&st, volume, nx, &v_emlrtRTEI);
  volume_data = volume->data;
  if (loop_ub_tmp != n) {
    emlrtIntegerCheckR2012b(imax, &c_emlrtDCI, &st);
  }
  for (n = 0; n < loop_ub_tmp; n++) {
    volume_data[n] = 0.0;
  }
  arc_length = 0.0;
  n = loop_ub_tmp - 1;
  emlrtForLoopVectorCheckR2021a(1.0, 1.0, imax - 1.0, mxDOUBLE_CLASS,
                                (int32_T)imax - 1, &c_emlrtRTEI, &st);
  for (i = 0; i < n; i++) {
    b_st.site = &i_emlrtRSI;
    if ((int32_T)((uint32_T)i + 1U) > x->size[0]) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)i + 1U), 1, x->size[0],
                                    &hc_emlrtBCI, &b_st);
    }
    if ((int32_T)((uint32_T)i + 2U) > x->size[0]) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)i + 2U), 1, x->size[0],
                                    &ic_emlrtBCI, &b_st);
    }
    dl = x_data[i + 1] - x_data[i];
    b_st.site = &i_emlrtRSI;
    if ((int32_T)((uint32_T)i + 1U) > y->size[0]) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)i + 1U), 1, y->size[0],
                                    &dc_emlrtBCI, &b_st);
    }
    if ((int32_T)((uint32_T)i + 2U) > y->size[0]) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)i + 2U), 1, y->size[0],
                                    &ec_emlrtBCI, &b_st);
    }
    vblend = y_data[i + 1] - y_data[i];
    dl = dl * dl + vblend * vblend;
    b_st.site = &i_emlrtRSI;
    if (dl < 0.0) {
      emlrtErrorWithMessageIdR2018a(
          &b_st, &u_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
          "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
    }
    dl = muDoubleScalarSqrt(dl);
    arc_length += dl;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(&st);
    }
  }
  n = loop_ub_tmp - 2;
  emlrtForLoopVectorCheckR2021a(2.0, 1.0, imax - 1.0, mxDOUBLE_CLASS,
                                (int32_T)imax - 2, &d_emlrtRTEI, &st);
  for (i = 0; i < n; i++) {
    if (i + 3 > rj->size[0]) {
      emlrtDynamicBoundsCheckR2012b(i + 3, 1, rj->size[0], &jc_emlrtBCI, &st);
    }
    if (i + 3 > rjm1->size[0]) {
      emlrtDynamicBoundsCheckR2012b(i + 3, 1, rjm1->size[0], &kc_emlrtBCI, &st);
    }
    if (i + 1 > rj->size[0]) {
      emlrtDynamicBoundsCheckR2012b(i + 1, 1, rj->size[0], &lc_emlrtBCI, &st);
    }
    if (i + 1 > rjm1->size[0]) {
      emlrtDynamicBoundsCheckR2012b(i + 1, 1, rjm1->size[0], &mc_emlrtBCI, &st);
    }
    b_st.site = &j_emlrtRSI;
    if ((int32_T)((uint32_T)i + 1U) > x->size[0]) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)i + 1U), 1, x->size[0],
                                    &bc_emlrtBCI, &b_st);
    }
    if ((int32_T)((uint32_T)i + 3U) > x->size[0]) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)i + 3U), 1, x->size[0],
                                    &cc_emlrtBCI, &b_st);
    }
    dl = (x_data[i + 2] - x_data[i]) / 2.0;
    b_st.site = &j_emlrtRSI;
    if ((int32_T)((uint32_T)i + 1U) > y->size[0]) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)i + 1U), 1, y->size[0],
                                    &yb_emlrtBCI, &b_st);
    }
    if ((int32_T)((uint32_T)i + 3U) > y->size[0]) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)i + 3U), 1, y->size[0],
                                    &ac_emlrtBCI, &b_st);
    }
    vblend = (y_data[i + 2] - y_data[i]) / 2.0;
    dl = dl * dl + vblend * vblend;
    b_st.site = &j_emlrtRSI;
    if (dl < 0.0) {
      emlrtErrorWithMessageIdR2018a(
          &b_st, &u_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
          "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
    }
    dl = muDoubleScalarSqrt(dl);
    if ((int32_T)((uint32_T)i + 2U) > volume->size[0]) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)i + 2U), 1,
                                    volume->size[0], &ib_emlrtBCI, &st);
    }
    if ((int32_T)((uint32_T)i + 2U) > scale->size[0]) {
      emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)i + 2U), 1,
                                    scale->size[0], &jb_emlrtBCI, &st);
    }
    vblend = scale_data[i + 1];
    volume_data[i + 1] =
        (((rj_data[i + 2] - rjm1_data[i + 2]) + rj_data[i]) - rjm1_data[i]) /
        2.0 * (vblend * dl + (1.0 - vblend) * arc_length / (jmax - 1.0));
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(&st);
    }
  }
  if (volume->size[0] < 1) {
    emlrtDynamicBoundsCheckR2012b(1, 1, volume->size[0], &oc_emlrtBCI, &st);
  }
  if (volume->size[0] < 2) {
    emlrtDynamicBoundsCheckR2012b(2, 1, volume->size[0], &nc_emlrtBCI, &st);
  }
  volume_data[0] = volume_data[1];
  if ((imax < 1.0) || (loop_ub_tmp > volume->size[0])) {
    emlrtDynamicBoundsCheckR2012b((int32_T)imax, 1, volume->size[0],
                                  &fc_emlrtBCI, &st);
  }
  if ((imax - 1.0 < 1.0) || (loop_ub_tmp - 1 > volume->size[0])) {
    emlrtDynamicBoundsCheckR2012b((int32_T)imax - 1, 1, volume->size[0],
                                  &gc_emlrtBCI, &st);
  }
  volume_data[loop_ub_tmp - 1] = volume_data[loop_ub_tmp - 2];
  /*  Calculate the xi metrics and estimate the eta metrics */
  emxInit_real_T(sp, &xxi, 1, &bb_emlrtRTEI);
  emxInit_real_T(sp, &yxi, 1, &bb_emlrtRTEI);
  emxInit_real_T(sp, &xeta, 1, &bb_emlrtRTEI);
  emxInit_real_T(sp, &yeta, 1, &bb_emlrtRTEI);
  st.site = &b_emlrtRSI;
  grid_metrics(&st, imax, x, y, volume, xxi, yxi, xeta, yeta);
  x_update_data = yeta->data;
  xeta_data = xeta->data;
  yxi_data = yxi->data;
  xxi_data = xxi->data;
  /*  Set up the matrix problem */
  st.site = &c_emlrtRSI;
  /* CREATE_RHS Creates the standard rhs matrix for volume marching plus */
  /*  See Steger & Chaussee eq's 10, 11, and 12 */
  emxInit_real_T(&st, &rhs, 1, &w_emlrtRTEI);
  u = 2U * (uint32_T)imax;
  nx = (int32_T)u;
  if ((real_T)u != (int32_T)u) {
    emlrtIntegerCheckR2012b(u, &emlrtDCI, &st);
  }
  n = rhs->size[0];
  rhs->size[0] = (int32_T)u;
  emxEnsureCapacity_real_T(&st, rhs, n, &w_emlrtRTEI);
  rhs_data = rhs->data;
  for (n = 0; n < nx; n++) {
    rhs_data[n] = 0.0;
  }
  /*  Compute RHS for Interior Points */
  for (i = 0; i < loop_ub_tmp; i++) {
    real_T dissipation_x;
    /*  Explicit fourth order dissipation */
    if (i + 1 == 1) {
      if (mu->size[0] < 1) {
        emlrtDynamicBoundsCheckR2012b(1, 1, mu->size[0], &pc_emlrtBCI, &st);
      }
      if (x->size[0] < 1) {
        emlrtDynamicBoundsCheckR2012b(1, 1, x->size[0], &qc_emlrtBCI, &st);
      }
      if (x->size[0] < 2) {
        emlrtDynamicBoundsCheckR2012b(2, 1, x->size[0], &rc_emlrtBCI, &st);
      }
      if (x->size[0] < 3) {
        emlrtDynamicBoundsCheckR2012b(3, 1, x->size[0], &sc_emlrtBCI, &st);
      }
      if (x->size[0] < 4) {
        emlrtDynamicBoundsCheckR2012b(4, 1, x->size[0], &tc_emlrtBCI, &st);
      }
      if (x->size[0] < 5) {
        emlrtDynamicBoundsCheckR2012b(5, 1, x->size[0], &uc_emlrtBCI, &st);
      }
      if (x->size[0] < 6) {
        emlrtDynamicBoundsCheckR2012b(6, 1, x->size[0], &vc_emlrtBCI, &st);
      }
      dissipation_x =
          -mu_data[0] *
          (((((3.0 * x_data[0] - 14.0 * x_data[1]) + 26.0 * x_data[2]) -
             24.0 * x_data[3]) +
            11.0 * x_data[4]) -
           2.0 * x_data[5]);
      if (mu->size[0] < 1) {
        emlrtDynamicBoundsCheckR2012b(1, 1, mu->size[0], &wc_emlrtBCI, &st);
      }
      if (y->size[0] < 1) {
        emlrtDynamicBoundsCheckR2012b(1, 1, y->size[0], &xc_emlrtBCI, &st);
      }
      if (y->size[0] < 2) {
        emlrtDynamicBoundsCheckR2012b(2, 1, y->size[0], &yc_emlrtBCI, &st);
      }
      if (y->size[0] < 3) {
        emlrtDynamicBoundsCheckR2012b(3, 1, y->size[0], &ad_emlrtBCI, &st);
      }
      if (y->size[0] < 4) {
        emlrtDynamicBoundsCheckR2012b(4, 1, y->size[0], &bd_emlrtBCI, &st);
      }
      if (y->size[0] < 5) {
        emlrtDynamicBoundsCheckR2012b(5, 1, y->size[0], &cd_emlrtBCI, &st);
      }
      if (y->size[0] < 6) {
        emlrtDynamicBoundsCheckR2012b(6, 1, y->size[0], &dd_emlrtBCI, &st);
      }
      dl = -mu_data[0] *
           (((((3.0 * y_data[0] - 14.0 * y_data[1]) + 26.0 * y_data[2]) -
              24.0 * y_data[3]) +
             11.0 * y_data[4]) -
            2.0 * y_data[5]);
    } else if (i + 1 == loop_ub_tmp) {
      if ((i - 4 < 1) || (i - 4 > x->size[0])) {
        emlrtDynamicBoundsCheckR2012b(i - 4, 1, x->size[0], &rb_emlrtBCI, &st);
      }
      if ((i - 3 < 1) || (i - 3 > x->size[0])) {
        emlrtDynamicBoundsCheckR2012b(i - 3, 1, x->size[0], &sb_emlrtBCI, &st);
      }
      if ((i - 2 < 1) || (i - 2 > x->size[0])) {
        emlrtDynamicBoundsCheckR2012b(i - 2, 1, x->size[0], &tb_emlrtBCI, &st);
      }
      if ((i - 1 < 1) || (i - 1 > x->size[0])) {
        emlrtDynamicBoundsCheckR2012b(i - 1, 1, x->size[0], &ub_emlrtBCI, &st);
      }
      if ((i < 1) || (i > x->size[0])) {
        emlrtDynamicBoundsCheckR2012b(i, 1, x->size[0], &vb_emlrtBCI, &st);
      }
      if (i + 1 > x->size[0]) {
        emlrtDynamicBoundsCheckR2012b(i + 1, 1, x->size[0], &wb_emlrtBCI, &st);
      }
      if (i + 1 > mu->size[0]) {
        emlrtDynamicBoundsCheckR2012b(i + 1, 1, mu->size[0], &xb_emlrtBCI, &st);
      }
      dl = -mu_data[i];
      dissipation_x =
          dl *
          (((((3.0 * x_data[i] - 14.0 * x_data[i - 1]) + 26.0 * x_data[i - 2]) -
             24.0 * x_data[i - 3]) +
            11.0 * x_data[i - 4]) -
           2.0 * x_data[i - 5]);
      if ((i - 4 < 1) || (i - 4 > y->size[0])) {
        emlrtDynamicBoundsCheckR2012b(i - 4, 1, y->size[0], &kb_emlrtBCI, &st);
      }
      if ((i - 3 < 1) || (i - 3 > y->size[0])) {
        emlrtDynamicBoundsCheckR2012b(i - 3, 1, y->size[0], &lb_emlrtBCI, &st);
      }
      if ((i - 2 < 1) || (i - 2 > y->size[0])) {
        emlrtDynamicBoundsCheckR2012b(i - 2, 1, y->size[0], &mb_emlrtBCI, &st);
      }
      if ((i - 1 < 1) || (i - 1 > y->size[0])) {
        emlrtDynamicBoundsCheckR2012b(i - 1, 1, y->size[0], &nb_emlrtBCI, &st);
      }
      if (i > y->size[0]) {
        emlrtDynamicBoundsCheckR2012b(i, 1, y->size[0], &ob_emlrtBCI, &st);
      }
      if (i + 1 > y->size[0]) {
        emlrtDynamicBoundsCheckR2012b(i + 1, 1, y->size[0], &pb_emlrtBCI, &st);
      }
      if (i + 1 > mu->size[0]) {
        emlrtDynamicBoundsCheckR2012b(i + 1, 1, mu->size[0], &qb_emlrtBCI, &st);
      }
      dl *=
          ((((3.0 * y_data[i] - 14.0 * y_data[i - 1]) + 26.0 * y_data[i - 2]) -
            24.0 * y_data[i - 3]) +
           11.0 * y_data[i - 4]) -
          2.0 * y_data[i - 5];
    } else if (i + 1 == 2) {
      if (mu->size[0] < 2) {
        emlrtDynamicBoundsCheckR2012b(2, 1, mu->size[0], &ed_emlrtBCI, &st);
      }
      if (x->size[0] < 1) {
        emlrtDynamicBoundsCheckR2012b(1, 1, x->size[0], &fd_emlrtBCI, &st);
      }
      if (x->size[0] < 2) {
        emlrtDynamicBoundsCheckR2012b(2, 1, x->size[0], &gd_emlrtBCI, &st);
      }
      if (x->size[0] < 3) {
        emlrtDynamicBoundsCheckR2012b(3, 1, x->size[0], &hd_emlrtBCI, &st);
      }
      if (x->size[0] < 4) {
        emlrtDynamicBoundsCheckR2012b(4, 1, x->size[0], &id_emlrtBCI, &st);
      }
      if (x->size[0] < 5) {
        emlrtDynamicBoundsCheckR2012b(5, 1, x->size[0], &jd_emlrtBCI, &st);
      }
      if (x->size[0] < 6) {
        emlrtDynamicBoundsCheckR2012b(6, 1, x->size[0], &kd_emlrtBCI, &st);
      }
      dissipation_x =
          -mu_data[1] *
          (((((2.0 * x_data[0] - 9.0 * x_data[1]) + 16.0 * x_data[2]) -
             14.0 * x_data[3]) +
            6.0 * x_data[4]) -
           x_data[5]);
      if (mu->size[0] < 2) {
        emlrtDynamicBoundsCheckR2012b(2, 1, mu->size[0], &ld_emlrtBCI, &st);
      }
      if (y->size[0] < 1) {
        emlrtDynamicBoundsCheckR2012b(1, 1, y->size[0], &md_emlrtBCI, &st);
      }
      if (y->size[0] < 2) {
        emlrtDynamicBoundsCheckR2012b(2, 1, y->size[0], &nd_emlrtBCI, &st);
      }
      if (y->size[0] < 3) {
        emlrtDynamicBoundsCheckR2012b(3, 1, y->size[0], &od_emlrtBCI, &st);
      }
      if (y->size[0] < 4) {
        emlrtDynamicBoundsCheckR2012b(4, 1, y->size[0], &pd_emlrtBCI, &st);
      }
      if (y->size[0] < 5) {
        emlrtDynamicBoundsCheckR2012b(5, 1, y->size[0], &qd_emlrtBCI, &st);
      }
      if (y->size[0] < 6) {
        emlrtDynamicBoundsCheckR2012b(6, 1, y->size[0], &rd_emlrtBCI, &st);
      }
      dl = -mu_data[1] *
           (((((2.0 * y_data[0] - 9.0 * y_data[1]) + 16.0 * y_data[2]) -
              14.0 * y_data[3]) +
             6.0 * y_data[4]) -
            y_data[5]);
    } else if (i + 1 == loop_ub_tmp - 1) {
      if ((i - 3 < 1) || (i - 3 > x->size[0])) {
        emlrtDynamicBoundsCheckR2012b(i - 3, 1, x->size[0], &t_emlrtBCI, &st);
      }
      if ((i - 2 < 1) || (i - 2 > x->size[0])) {
        emlrtDynamicBoundsCheckR2012b(i - 2, 1, x->size[0], &u_emlrtBCI, &st);
      }
      if ((i - 1 < 1) || (i - 1 > x->size[0])) {
        emlrtDynamicBoundsCheckR2012b(i - 1, 1, x->size[0], &v_emlrtBCI, &st);
      }
      if ((i < 1) || (i > x->size[0])) {
        emlrtDynamicBoundsCheckR2012b(i, 1, x->size[0], &w_emlrtBCI, &st);
      }
      if (i + 1 > x->size[0]) {
        emlrtDynamicBoundsCheckR2012b(i + 1, 1, x->size[0], &x_emlrtBCI, &st);
      }
      if ((int32_T)((uint32_T)i + 2U) > x->size[0]) {
        emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)i + 2U), 1,
                                      x->size[0], &y_emlrtBCI, &st);
      }
      if (i + 1 > mu->size[0]) {
        emlrtDynamicBoundsCheckR2012b(i + 1, 1, mu->size[0], &ab_emlrtBCI, &st);
      }
      dl = -mu_data[i];
      dissipation_x =
          dl *
          (((((2.0 * x_data[i + 1] - 9.0 * x_data[i]) + 16.0 * x_data[i - 1]) -
             14.0 * x_data[i - 2]) +
            6.0 * x_data[i - 3]) -
           x_data[i - 4]);
      if ((i - 3 < 1) || (i - 3 > y->size[0])) {
        emlrtDynamicBoundsCheckR2012b(i - 3, 1, y->size[0], &f_emlrtBCI, &st);
      }
      if ((i - 2 < 1) || (i - 2 > y->size[0])) {
        emlrtDynamicBoundsCheckR2012b(i - 2, 1, y->size[0], &g_emlrtBCI, &st);
      }
      if ((i - 1 < 1) || (i - 1 > y->size[0])) {
        emlrtDynamicBoundsCheckR2012b(i - 1, 1, y->size[0], &h_emlrtBCI, &st);
      }
      if (i > y->size[0]) {
        emlrtDynamicBoundsCheckR2012b(i, 1, y->size[0], &i_emlrtBCI, &st);
      }
      if (i + 1 > y->size[0]) {
        emlrtDynamicBoundsCheckR2012b(i + 1, 1, y->size[0], &j_emlrtBCI, &st);
      }
      if ((int32_T)((uint32_T)i + 2U) > y->size[0]) {
        emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)i + 2U), 1,
                                      y->size[0], &k_emlrtBCI, &st);
      }
      if (i + 1 > mu->size[0]) {
        emlrtDynamicBoundsCheckR2012b(i + 1, 1, mu->size[0], &l_emlrtBCI, &st);
      }
      dl *= ((((2.0 * y_data[i + 1] - 9.0 * y_data[i]) + 16.0 * y_data[i - 1]) -
              14.0 * y_data[i - 2]) +
             6.0 * y_data[i - 3]) -
            y_data[i - 4];
    } else {
      if (((int32_T)((uint32_T)i + 3U) < 1) ||
          ((int32_T)((uint32_T)i + 3U) > x->size[0])) {
        emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)i + 3U), 1,
                                      x->size[0], &bb_emlrtBCI, &st);
      }
      if (((int32_T)((uint32_T)i + 2U) < 1) ||
          ((int32_T)((uint32_T)i + 2U) > x->size[0])) {
        emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)i + 2U), 1,
                                      x->size[0], &cb_emlrtBCI, &st);
      }
      if (i + 1 > x->size[0]) {
        emlrtDynamicBoundsCheckR2012b(i + 1, 1, x->size[0], &db_emlrtBCI, &st);
      }
      if ((i < 1) || (i > x->size[0])) {
        emlrtDynamicBoundsCheckR2012b(i, 1, x->size[0], &eb_emlrtBCI, &st);
      }
      if ((i - 1 < 1) || (i - 1 > x->size[0])) {
        emlrtDynamicBoundsCheckR2012b(i - 1, 1, x->size[0], &fb_emlrtBCI, &st);
      }
      if (i + 1 > mu->size[0]) {
        emlrtDynamicBoundsCheckR2012b(i + 1, 1, mu->size[0], &gb_emlrtBCI, &st);
      }
      dl = -mu_data[i];
      dissipation_x =
          dl * ((((x_data[i - 2] - 4.0 * x_data[i - 1]) + 6.0 * x_data[i]) -
                 4.0 * x_data[i + 1]) +
                x_data[i + 2]);
      if ((int32_T)((uint32_T)i + 3U) > y->size[0]) {
        emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)i + 3U), 1,
                                      y->size[0], &m_emlrtBCI, &st);
      }
      if ((int32_T)((uint32_T)i + 2U) > y->size[0]) {
        emlrtDynamicBoundsCheckR2012b((int32_T)((uint32_T)i + 2U), 1,
                                      y->size[0], &n_emlrtBCI, &st);
      }
      if (i + 1 > y->size[0]) {
        emlrtDynamicBoundsCheckR2012b(i + 1, 1, y->size[0], &o_emlrtBCI, &st);
      }
      if (i > y->size[0]) {
        emlrtDynamicBoundsCheckR2012b(i, 1, y->size[0], &p_emlrtBCI, &st);
      }
      if ((i - 1 < 1) || (i - 1 > y->size[0])) {
        emlrtDynamicBoundsCheckR2012b(i - 1, 1, y->size[0], &q_emlrtBCI, &st);
      }
      if (i + 1 > mu->size[0]) {
        emlrtDynamicBoundsCheckR2012b(i + 1, 1, mu->size[0], &r_emlrtBCI, &st);
      }
      dl *= (((y_data[i - 2] - 4.0 * y_data[i - 1]) + 6.0 * y_data[i]) -
             4.0 * y_data[i + 1]) +
            y_data[i + 2];
    }
    if (i + 1 > xxi->size[0]) {
      emlrtDynamicBoundsCheckR2012b(i + 1, 1, xxi->size[0], &sd_emlrtBCI, &st);
    }
    if (i + 1 > yxi->size[0]) {
      emlrtDynamicBoundsCheckR2012b(i + 1, 1, yxi->size[0], &td_emlrtBCI, &st);
    }
    if (i + 1 > xeta->size[0]) {
      emlrtDynamicBoundsCheckR2012b(i + 1, 1, xeta->size[0], &ud_emlrtBCI, &st);
    }
    if (i + 1 > yeta->size[0]) {
      emlrtDynamicBoundsCheckR2012b(i + 1, 1, yeta->size[0], &vd_emlrtBCI, &st);
    }
    arc_length = 1.0 / (xxi_data[i] * xxi_data[i] + yxi_data[i] * yxi_data[i]);
    if (i + 1 > scale->size[0]) {
      emlrtDynamicBoundsCheckR2012b(i + 1, 1, scale->size[0], &c_emlrtBCI, &st);
    }
    if (i + 1 > volume->size[0]) {
      emlrtDynamicBoundsCheckR2012b(i + 1, 1, volume->size[0], &d_emlrtBCI,
                                    &st);
    }
    if (i + 1 > scale->size[0]) {
      emlrtDynamicBoundsCheckR2012b(i + 1, 1, scale->size[0], &e_emlrtBCI, &st);
    }
    vblend = scale_data[i];
    vblend = vblend * volume_data[i] +
             (1.0 - vblend) *
                 (xxi_data[i] * x_update_data[i] - yxi_data[i] * xeta_data[i]);
    u = (uint32_T)i << 1;
    if (((int32_T)(u + 1U) < 1) || ((int32_T)(u + 1U) > rhs->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(u + 1U), 1, rhs->size[0],
                                    &b_emlrtBCI, &st);
    }
    rhs_data[(int32_T)u] = -yxi_data[i] * arc_length * vblend + dissipation_x;
    if (((int32_T)(u + 2U) < 1) || ((int32_T)(u + 2U) > rhs->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(u + 2U), 1, rhs->size[0],
                                    &emlrtBCI, &st);
    }
    rhs_data[(int32_T)u + 1] = xxi_data[i] * arc_length * vblend + dl;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(&st);
    }
  }
  /*  Compute BC */
  if (rhs->size[0] < 1) {
    emlrtDynamicBoundsCheckR2012b(1, 1, rhs->size[0], &wd_emlrtBCI, &st);
  }
  rhs_data[0] = 0.0;
  if (rhs->size[0] < 2) {
    emlrtDynamicBoundsCheckR2012b(2, 1, rhs->size[0], &xd_emlrtBCI, &st);
  }
  rhs_data[1] = 0.0;
  u = (uint32_T)(loop_ub_tmp - 1) << 1;
  if (((int32_T)(u + 1U) < 1) || ((int32_T)(u + 1U) > rhs->size[0])) {
    emlrtDynamicBoundsCheckR2012b((int32_T)(u + 1U), 1, rhs->size[0],
                                  &hb_emlrtBCI, &st);
  }
  rhs_data[(int32_T)u] = 0.0;
  if (((int32_T)(u + 2U) < 1) || ((int32_T)(u + 2U) > rhs->size[0])) {
    emlrtDynamicBoundsCheckR2012b((int32_T)(u + 2U), 1, rhs->size[0],
                                  &s_emlrtBCI, &st);
  }
  rhs_data[(int32_T)u + 1] = 0.0;
  emxInit_real_T(sp, &lhs_d, 1, &ab_emlrtRTEI);
  emxInit_int32_T(sp, &lhs_colidx, &ab_emlrtRTEI);
  emxInit_int32_T(sp, &lhs_rowidx, &ab_emlrtRTEI);
  st.site = &d_emlrtRSI;
  nx = create_lhs(&st, imax, xxi, yxi, xeta, yeta, alpham, muim, lhs_d,
                  lhs_colidx, lhs_rowidx, &i);
  emxFree_real_T(sp, &yeta);
  emxFree_real_T(sp, &xeta);
  emxFree_real_T(sp, &yxi);
  emxFree_real_T(sp, &xxi);
  /*  Solve and update */
  st.site = &e_emlrtRSI;
  sparse_mldivide(&st, lhs_d, lhs_colidx, lhs_rowidx, nx, i, rhs, volume);
  volume_data = volume->data;
  emxFree_int32_T(sp, &lhs_rowidx);
  emxFree_int32_T(sp, &lhs_colidx);
  emxFree_real_T(sp, &lhs_d);
  emxFree_real_T(sp, &rhs);
  dv[0] = 2.0;
  dv[1] = imax;
  st.site = &f_emlrtRSI;
  nx = volume->size[0];
  b_st.site = &dc_emlrtRSI;
  assertValidSizeArg(&b_st, dv);
  n = volume->size[0];
  if (volume->size[0] < 1) {
    n = 1;
  }
  nx = muIntScalarMax_sint32(nx, n);
  if (nx < 2) {
    emlrtErrorWithMessageIdR2018a(&st, &emlrtRTEI,
                                  "Coder:toolbox:reshape_emptyReshapeLimit",
                                  "Coder:toolbox:reshape_emptyReshapeLimit", 0);
  }
  if (loop_ub_tmp > nx) {
    emlrtErrorWithMessageIdR2018a(&st, &emlrtRTEI,
                                  "Coder:toolbox:reshape_emptyReshapeLimit",
                                  "Coder:toolbox:reshape_emptyReshapeLimit", 0);
  }
  if (loop_ub_tmp << 1 != volume->size[0]) {
    emlrtErrorWithMessageIdR2018a(
        &st, &b_emlrtRTEI, "Coder:MATLAB:getReshapeDims_notSameNumel",
        "Coder:MATLAB:getReshapeDims_notSameNumel", 0);
  }
  dv[0] = imax;
  dv[1] = 1.0;
  st.site = &g_emlrtRSI;
  b_st.site = &dc_emlrtRSI;
  assertValidSizeArg(&b_st, dv);
  n = 1;
  if (imax > 1.0) {
    n = (int32_T)imax;
  }
  if (loop_ub_tmp > muIntScalarMax_sint32(loop_ub_tmp, n)) {
    emlrtErrorWithMessageIdR2018a(&st, &emlrtRTEI,
                                  "Coder:toolbox:reshape_emptyReshapeLimit",
                                  "Coder:toolbox:reshape_emptyReshapeLimit", 0);
  }
  if ((x->size[0] != loop_ub_tmp) && ((x->size[0] != 1) && (imax != 1.0))) {
    emlrtDimSizeImpxCheckR2021b(x->size[0], (int32_T)imax, &emlrtECI,
                                (emlrtConstCTX)sp);
  }
  if (x->size[0] == loop_ub_tmp) {
    n = x_update->size[0];
    x_update->size[0] = x->size[0];
    emxEnsureCapacity_real_T(sp, x_update, n, &x_emlrtRTEI);
    x_update_data = x_update->data;
    nx = x->size[0];
    for (n = 0; n < nx; n++) {
      x_update_data[n] = x_data[n] + volume_data[2 * n];
    }
  } else {
    st.site = &g_emlrtRSI;
    b_binary_expand_op(&st, x_update, x, volume, dv);
  }
  st.site = &h_emlrtRSI;
  i = (int32_T)imax;
  b_st.site = &dc_emlrtRSI;
  assertValidSizeArg(&b_st, dv);
  n = 1;
  if (imax > 1.0) {
    n = (int32_T)imax;
  }
  if (loop_ub_tmp > muIntScalarMax_sint32(i, n)) {
    emlrtErrorWithMessageIdR2018a(&st, &emlrtRTEI,
                                  "Coder:toolbox:reshape_emptyReshapeLimit",
                                  "Coder:toolbox:reshape_emptyReshapeLimit", 0);
  }
  if (imax != imax) {
    emlrtErrorWithMessageIdR2018a(
        &st, &b_emlrtRTEI, "Coder:MATLAB:getReshapeDims_notSameNumel",
        "Coder:MATLAB:getReshapeDims_notSameNumel", 0);
  }
  n = (int32_T)imax;
  if ((y->size[0] != n) && ((y->size[0] != 1) && (imax != 1.0))) {
    emlrtDimSizeImpxCheckR2021b(y->size[0], (int32_T)imax, &b_emlrtECI,
                                (emlrtConstCTX)sp);
  }
  if (y->size[0] == n) {
    n = y_update->size[0];
    y_update->size[0] = y->size[0];
    emxEnsureCapacity_real_T(sp, y_update, n, &y_emlrtRTEI);
    x_update_data = y_update->data;
    nx = y->size[0];
    for (n = 0; n < nx; n++) {
      x_update_data[n] = y_data[n] + volume_data[2 * n + 1];
    }
  } else {
    st.site = &h_emlrtRSI;
    binary_expand_op(&st, y_update, y, volume, dv);
  }
  emxFree_real_T(sp, &volume);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

/* End of code generation (march_grid_fast.c) */
