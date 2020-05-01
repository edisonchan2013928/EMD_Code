#pragma once
#ifndef UB_H_H
#define UB_H_H

#include "init.h"

void initHilbert_File(char*hilbertDimVector_FileName, emd_struct& my_emd);
double UB_hilbert(double*q, double*p, emd_struct& my_emd);
void UB_Hilbert_Pair_Value_Calculation(emd_struct& my_emd);

#endif