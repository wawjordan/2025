/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_march_grid_fast_mex.c
 *
 * Code generation for function '_coder_march_grid_fast_mex'
 *
 */

/* Include files */
#include "_coder_march_grid_fast_mex.h"
#include "_coder_march_grid_fast_api.h"
#include "march_grid_fast_data.h"
#include "march_grid_fast_initialize.h"
#include "march_grid_fast_terminate.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void march_grid_fast_mexFunction(int32_T nlhs, mxArray *plhs[2], int32_T nrhs,
                                 const mxArray *prhs[10])
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  const mxArray *outputs[2];
  int32_T i;
  st.tls = emlrtRootTLSGlobal;
  /* Check for proper number of arguments. */
  if (nrhs != 10) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 10, 4,
                        15, "march_grid_fast");
  }
  if (nlhs > 2) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 15,
                        "march_grid_fast");
  }
  /* Call the function. */
  march_grid_fast_api(prhs, nlhs, outputs);
  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    i = 1;
  } else {
    i = nlhs;
  }
  emlrtReturnArrays(i, &plhs[0], &outputs[0]);
}

void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs,
                 const mxArray *prhs[])
{
  mexAtExit(&march_grid_fast_atexit);
  /* Module initialization. */
  march_grid_fast_initialize();
  /* Dispatch the entry-point. */
  march_grid_fast_mexFunction(nlhs, plhs, nrhs, prhs);
  /* Module termination. */
  march_grid_fast_terminate();
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLSR2022a(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1,
                           NULL, "windows-1252", true);
  return emlrtRootTLSGlobal;
}

/* End of code generation (_coder_march_grid_fast_mex.c) */
