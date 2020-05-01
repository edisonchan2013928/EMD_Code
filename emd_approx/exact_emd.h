#pragma once
#ifndef EXACT_EMD_H
#define EXACT_EMD_H
#include "init.h"
#include "emd.h"

void EMD_PairCalculation(emd_struct& my_emd);
void costMatrix_Reduction(double*q, double*p, int alpha, int beta, double**costMatrix, int dim);
void IgnoreZero(double*q, double*p, vector<double>& q_bar, vector<double>& p_bar, double**costMatrix, int dim);
double emdComputation(vector<double>& q_bar, vector<double>& p_bar);
void clearCostMatrix(int q_DIM_Num);

#endif