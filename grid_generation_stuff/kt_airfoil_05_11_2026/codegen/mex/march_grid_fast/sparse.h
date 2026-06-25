/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sparse.h
 *
 * Code generation for function 'sparse'
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
void sparse_mldivide(const emlrtStack *sp, const emxArray_real_T *A_d,
                     const emxArray_int32_T *A_colidx,
                     const emxArray_int32_T *A_rowidx, int32_T A_m, int32_T A_n,
                     const emxArray_real_T *b, emxArray_real_T *y);

/* End of code generation (sparse.h) */
