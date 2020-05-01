#include "UB_G.h"

const double epsilon = 0.00001;

bool compare(index_Cost lhs, index_Cost rhs)
{
	return lhs.costValue < rhs.costValue;
}

void list_PreProcess(double**costMatrix, vector<index_Cost>& sortCostList, int dim)
{
	//clear for safety
	sortCostList.clear();

	index_Cost temp;
	for (int i = 0; i < dim; i++)
	{
		temp.qIndex = i;
		for (int j = 0; j < dim; j++)
		{
			temp.pIndex = j;
			temp.costValue = costMatrix[i][j];
			sortCostList.push_back(temp);
		}
	}
	//Sort it in increasing order by value
	sort(sortCostList.begin(), sortCostList.end(), compare);
}

inline void sparseProcessing(vector<index_Cost>& sortCostList, double**costMatrix, double*qVector, double*pVector, int dim, int& dim1, int& dim2)
{
	vector<int> qRedIndexList;
	vector<int> pRedIndexList;
	index_Cost temp;

	sortCostList.clear();
	dim1 = 0, dim2 = 0;

	for (int i = 0; i < dim; i++)
	{
		if (qVector[i] > epsilon)
		{
			qRedIndexList.push_back(i);
			dim1++;
		}

		if (pVector[i] > epsilon)
		{
			pRedIndexList.push_back(i);
			dim2++;
		}
	}

	for (int ii = 0; ii < dim1; ii++)
	{
		temp.qIndex = qRedIndexList[ii];

		for (int jj = 0; jj < dim2; jj++)
		{
			temp.costValue = costMatrix[qRedIndexList[ii]][pRedIndexList[jj]];
			temp.pIndex = pRedIndexList[jj];

			sortCostList.push_back(temp);
		}
	}

	//Sort it in increasing order by value
	sort(sortCostList.begin(), sortCostList.end(), compare);

	qRedIndexList.clear();
	pRedIndexList.clear();
}

double UB_Greedy(double*qVector, double*pVector, emd_struct& my_emd)
{
	int dim_Square;

	for (int i = 0; i < my_emd.dim; i++)
	{
		my_emd.LVectorTemp[i] = qVector[i];
		my_emd.RVectorTemp[i] = pVector[i];
	}

	if (my_emd.isSparse == true)//Use Sparse Processing
	{
		int dim1;
		int dim2;
		sparseProcessing(my_emd.sortCostList, my_emd.costMatrix, qVector, pVector, my_emd.dim, dim1, dim2);
		dim_Square = dim1 * dim2;
	}
	else //No Sparse Processing
		dim_Square = my_emd.dim * my_emd.dim;

	double Cost = 0.0;
	double minFlow;

	for (int index = 0; index < dim_Square; index++)
	{
		minFlow = min(my_emd.LVectorTemp[my_emd.sortCostList[index].qIndex], my_emd.RVectorTemp[my_emd.sortCostList[index].pIndex]);
		Cost = Cost + my_emd.sortCostList[index].costValue*minFlow;
		my_emd.LVectorTemp[my_emd.sortCostList[index].qIndex] = my_emd.LVectorTemp[my_emd.sortCostList[index].qIndex] - minFlow;
		my_emd.RVectorTemp[my_emd.sortCostList[index].pIndex] = my_emd.RVectorTemp[my_emd.sortCostList[index].pIndex] - minFlow;
	}

	return Cost;
}

void UB_Greedy_Pair_Value_Calculation(emd_struct& my_emd)
{
	for (int p = 0; p < my_emd.pairNum; p++)
		my_emd.result_pairVector.push_back(UB_Greedy(my_emd.dataSetVector[my_emd.Pair_Array[p].L_Index], my_emd.dataSetVector[my_emd.Pair_Array[p].R_Index], my_emd));
}