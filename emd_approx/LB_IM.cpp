#include "LB_IM.h"

bool compFunctionCostAscend(component c1, component c2)
{
	return (c1.cost < c2.cost);
}

void init_component(emd_struct& my_emd)
{
	//init LVectorTemp and RVectorTemp
	my_emd.LVectorTemp = new double[my_emd.dim];
	my_emd.RVectorTemp = new double[my_emd.dim];

	vector<component> tempVectorComponent;
	component tempComponent;

	for (int i = 0; i < my_emd.dim; i++)
	{
		my_emd.SortCostVector.push_back(tempVectorComponent);
		for (int j = 0; j < my_emd.dim; j++)
		{
			tempComponent.rIndex = j;
			tempComponent.cost = my_emd.costMatrix[i][j];
			my_emd.SortCostVector[i].push_back(tempComponent);
		}
	}

	//Do the sorting in increasing order according to the cost value
	for (int i = 0; i < my_emd.dim; i++)
		sort(my_emd.SortCostVector[i].begin(), my_emd.SortCostVector[i].end(), compFunctionCostAscend);
}

double IM_LB(double*LVector, double*RVector, emd_struct& my_emd)
{
	double flow = 0.0;
	double IMValue = 0.0;

	for (int i = 0; i < my_emd.dim; i++)
	{
		my_emd.LVectorTemp[i] = LVector[i] - min(LVector[i], RVector[i]);
		my_emd.RVectorTemp[i] = RVector[i] - min(LVector[i], RVector[i]);
	}

	for (int i = 0; i < my_emd.dim; i++)
	{
		for (int j = 0; j < my_emd.dim; j++)
		{
			flow = min(my_emd.LVectorTemp[i], my_emd.RVectorTemp[my_emd.SortCostVector[i][j].rIndex]);
			IMValue = IMValue + flow * my_emd.SortCostVector[i][j].cost;
			my_emd.LVectorTemp[i] = my_emd.LVectorTemp[i] - flow;
			if (my_emd.LVectorTemp[i] < 0.0)
				continue;
		}
	}

	return IMValue;
}

void IM_Pair_Value_Calculation(emd_struct& my_emd)
{
	for (int p = 0; p < my_emd.pairNum; p++)
		my_emd.result_pairVector.push_back(IM_LB(my_emd.dataSetVector[my_emd.Pair_Array[p].L_Index], my_emd.dataSetVector[my_emd.Pair_Array[p].R_Index], my_emd));
}