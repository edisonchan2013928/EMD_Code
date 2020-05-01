#pragma once
#ifndef SINKHORN_H
#define SINKHORN_H

#include "math_operation.h"

void round(double**F_matrix, double*q, double*p, int dim);
void greenkhorn(double**A, double*q, double*p, int dim, double epsilon_dash);
void greenkhorn_Matlab(double**A, double*q, double*p, int dim, double epsilon_dash);
double approxOT(double**costMatrix, double*q, double*p, double abs_epsilon, int dim);
void optimal_transport_Sinkhorn(emd_struct& my_emd);

#endif