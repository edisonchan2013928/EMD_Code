#pragma once
#ifndef LB_DIM_RED_H
#define LB_DIM_RED_H

#include "init.h"
#include "exact_emd.h"

void randomGenerator_withoutReplacement(vector<int>& center, int dim_Red);
void globalMinimumClustering(double**costMatrix, vector<int>& center, vector< vector<int> >& center_dimList, int dim, int dim_Red);
double localMinimumClustering(vector<int>& center, vector< vector<int> >& center_dimList, int k, double**costMatrix);
void clearCenter_dimList(vector< vector<int> >& center_dimList, int dim_Red);
void init_dimRed(emd_struct& my_emd);
void kMetroid_Algorithm(emd_struct& my_emd);
void CostMatrixReduction_Metroid(emd_struct& my_emd);
void dataSetReduction(emd_struct& my_emd);
double LB_dimReduction(double*q_Red, double*p_Red, emd_struct& my_emd);
void LB_dimReduction_Pair_Value_Calculation(emd_struct& my_emd);

#endif