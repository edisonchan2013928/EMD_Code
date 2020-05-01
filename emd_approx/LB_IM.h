#pragma once
#ifndef LB_IM_H
#define LB_IM_H

#include "init.h"

void init_component(emd_struct& my_emd);
double IM_LB(double*LVector, double*RVector, emd_struct& my_emd);
void IM_Pair_Value_Calculation(emd_struct& my_emd);

#endif