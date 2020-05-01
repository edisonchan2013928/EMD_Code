#pragma once
#ifndef MATH_OPERATION_H
#define MATH_OPERATION_H

#include "init.h"

void row_sum(double**matrix, int dim, double*row_sumVec);
void col_sum(double**matrix, int dim, double*col_sumVec);
void outer_product(double*v1, double*v2, int dim, double**outerMatrix);
void matrix_divide_constant(double**matrix, int dim, double constant);
double L1_norm(double*vector, int dim);
double matrix_norm_L_1(double**matrix, int dim);
double matrix_norm_L_infty(double**matrix, int dim);
void exp_Matrix(double**C, double**A, double eta, int dim);
double matrix_inner_product(double**P, double**C, int dim);

#endif