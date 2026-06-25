/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * march_grid_fast.h
 *
 * Code generation for function 'march_grid_fast'
 *
 */

#pragma once

/* Include files */
#include "march_grid_fast_types.h"
#include "rtwtypes.h"
#include "emlrt.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
void march_grid_fast(const emlrtStack *sp, real_T imax, real_T jmax,
                     const emxArray_real_T *x, const emxArray_real_T *y,
                     const emxArray_real_T *mu, const emxArray_real_T *muim,
                     const emxArray_real_T *alpham,
                     const emxArray_real_T *scale, const emxArray_real_T *rj,
                     const emxArray_real_T *rjm1, emxArray_real_T *x_update,
                     emxArray_real_T *y_update);

/* End of code generation (march_grid_fast.h) */
