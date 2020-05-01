#include "LB_dimRed.h"

void randomGenerator_withoutReplacement(vector<int>& center, int dim_Red)
{
	srand(time(0));
	int gen;
	bool loop;

	gen = rand() % dim_Red;
	center.push_back(gen);

	for (int k = 1; k < dim_Red; k++)
	{
		loop = true;
		while (loop == true)
		{
			loop = false;
			gen = rand() % dim_Red;
			for (int g = 0; g < (int)center.size(); g++)
			{
				if (gen == center[g])
					loop = true;
			}
		}
		center.push_back(gen);
	}

}

void globalMinimumClustering(double**costMatrix, vector<int>& center, vector< vector<int> >& center_dimList, int dim, int dim_Red)
{
	double tempCost;
	int tempCenterIndex = -1;
	for (int d = 0; d < dim; d++)
	{
		tempCost = inf;
		tempCenterIndex = -1;
		for (int k = 0; k < dim_Red; k++)
		{
			if (costMatrix[d][center[k]] < tempCost)
			{
				tempCost = costMatrix[d][center[k]];
				tempCenterIndex = k;
			}
		}
		center_dimList[tempCenterIndex].push_back(d);
	}
}

double localMinimumClustering(vector<int>& center, vector< vector<int> >& center_dimList, int k, double**costMatrix)
{
	static vector<int> List;
	double minCost = inf;
	double tempCost;
	int centerIndex = -1;

	//List.push_back(center[k]);
	for (int i = 0; i < (int)center_dimList[k].size(); i++)
		List.push_back(center_dimList[k][i]);

	for (int c = 0; c < (int)List.size(); c++)
	{
		tempCost = 0;
		for (int l = 0; l < (int)List.size(); l++)
		{
			if (c != l)
				tempCost = tempCost + costMatrix[List[l]][List[c]];
		}
		if (tempCost < minCost)
		{
			minCost = tempCost;
			centerIndex = c;
		}
	}
	center[k] = List[centerIndex];
	//Do not need to update centerList, as it will be clear later 
	//or it will be reported when the cost does not change
	List.clear();

	return minCost;
}

void clearCenter_dimList(vector< vector<int> >& center_dimList, int dim_Red)
{
	for (int i = 0; i < dim_Red; i++)
		center_dimList[i].clear();
}

void init_dimRed(emd_struct& my_emd)
{
	vector<int>tempVector;
	for (int i = 0; i < my_emd.dim_Red; i++)
		my_emd.center_dimList.push_back(tempVector);

	my_emd.costMatrix_Red = new double*[my_emd.dim_Red];
	for (int d = 0; d < my_emd.dim_Red; d++)
		my_emd.costMatrix_Red[d] = new double[my_emd.dim_Red];

	my_emd.Red_dataSetVector = new double*[my_emd.dataNum];
	for (int i = 0; i < my_emd.dataNum; i++)
		my_emd.Red_dataSetVector[i] = new double[my_emd.dim_Red];

	for (int i = 0; i < my_emd.dataNum; i++)
		for (int d = 0; d < my_emd.dim_Red; d++)
			my_emd.Red_dataSetVector[i][d] = 0;
}

void kMetroid_Algorithm(emd_struct& my_emd)
{
	vector<int> center;
	bool converge = false;
	double TD = inf;
	double tempCost;

	randomGenerator_withoutReplacement(center, my_emd.dim_Red);

	while (converge != true)
	{
		tempCost = 0;
		globalMinimumClustering(my_emd.costMatrix, center, my_emd.center_dimList, my_emd.dim, my_emd.dim_Red);
		for (int k = 0; k < my_emd.dim_Red; k++)
			tempCost = tempCost + localMinimumClustering(center, my_emd.center_dimList, k, my_emd.costMatrix);

		if (tempCost < TD)
		{
			TD = tempCost;
			clearCenter_dimList(my_emd.center_dimList, my_emd.dim_Red);
		}
		else
			converge = true;
	}

}

void CostMatrixReduction_Metroid(emd_struct& my_emd)
{
	for (int m1 = 0; m1 < my_emd.dim_Red; m1++)
	{
		for (int m2 = 0; m2 < my_emd.dim_Red; m2++)
		{
			my_emd.costMatrix_Red[m1][m2] = inf;
			for (int i1 = 0; i1 < (int)my_emd.center_dimList[m1].size(); i1++)
			{
				for (int i2 = 0; i2 < (int)my_emd.center_dimList[m2].size(); i2++)
				{
					if (my_emd.costMatrix[my_emd.center_dimList[m1][i1]][my_emd.center_dimList[m2][i2]] < my_emd.costMatrix_Red[m1][m2])
						my_emd.costMatrix_Red[m1][m2] = my_emd.costMatrix[my_emd.center_dimList[m1][i1]][my_emd.center_dimList[m2][i2]];
				}
			}

			if (my_emd.costMatrix_Red[m1][m2] == inf)
				my_emd.costMatrix_Red[m1][m2] = 0;
		}
	}

}

void dataSetReduction(emd_struct& my_emd)
{
	double temp;
	for (int i = 0; i < my_emd.dataNum; i++)
	{
		for (int m = 0; m < my_emd.dim_Red; m++)
		{
			temp = 0;
			for (int ii = 0; ii < (int)my_emd.center_dimList[m].size(); ii++)
				temp = temp + my_emd.dataSetVector[i][my_emd.center_dimList[m][ii]];

			my_emd.Red_dataSetVector[i][m] = temp;
		}
	}
}

double LB_dimReduction(double*q_Red, double*p_Red, emd_struct& my_emd)
{
	double LB;
	IgnoreZero(q_Red, p_Red, my_emd.q_Red_bar, my_emd.p_Red_bar, my_emd.costMatrix_Red, my_emd.dim_Red);
	LB = emdComputation(my_emd.q_Red_bar, my_emd.p_Red_bar);

	return LB;
}

void LB_dimReduction_Pair_Value_Calculation(emd_struct& my_emd)
{
	for (int p = 0; p < my_emd.pairNum; p++)
		my_emd.result_pairVector.push_back(LB_dimReduction(my_emd.Red_dataSetVector[my_emd.Pair_Array[p].L_Index], my_emd.Red_dataSetVector[my_emd.Pair_Array[p].R_Index], my_emd));
}