#pragma once
#ifndef UB_G_H
#define UB_G_H

#include "init.h"

void list_PreProcess(double**costMatrix, vector<index_Cost>& sortCostList, int dim);
double UB_Greedy(double*qVector, double*pVector, emd_struct& my_emd);
void UB_Greedy_Pair_Value_Calculation(emd_struct& my_emd);

#endif