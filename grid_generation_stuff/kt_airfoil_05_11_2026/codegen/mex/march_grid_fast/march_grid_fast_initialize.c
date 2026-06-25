/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * march_grid_fast_initialize.c
 *
 * Code generation for function 'march_grid_fast_initialize'
 *
 */

/* Include files */
#include "march_grid_fast_initialize.h"
#include "_coder_march_grid_fast_mex.h"
#include "march_grid_fast_data.h"
#include "rt_nonfinite.h"

/* Function Declarations */
static void march_grid_fast_once(void);

/* Function Definitions */
static void march_grid_fast_once(void)
{
  mex_InitInfAndNan();
}

void march_grid_fast_initialize(void)
{
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
    march_grid_fast_once();
  }
}

/* End of code generation (march_grid_fast_initialize.c) */
