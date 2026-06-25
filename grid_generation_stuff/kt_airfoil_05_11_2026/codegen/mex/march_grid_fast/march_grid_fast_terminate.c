/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * march_grid_fast_terminate.c
 *
 * Code generation for function 'march_grid_fast_terminate'
 *
 */

/* Include files */
#include "march_grid_fast_terminate.h"
#include "_coder_march_grid_fast_mex.h"
#include "march_grid_fast_data.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void march_grid_fast_atexit(void)
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

void march_grid_fast_terminate(void)
{
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (march_grid_fast_terminate.c) */
