#include "UB_H.h"
#include "init.h"

void initHilbert_File(char*hilbertDimVector_FileName, emd_struct& my_emd)
{
	fstream hilbert_DimVector_File;
	hilbert_DimVector_File.open(hilbertDimVector_FileName);
	if (hilbert_DimVector_File.is_open() == false)
		cout << "Cannot open hilbert_DimVector_File!" << endl;

	my_emd.hilbert_DimVector = new int[my_emd.dim];
	for (int d = 0; d < my_emd.dim; d++)
		hilbert_DimVector_File >> my_emd.hilbert_DimVector[d];

	my_emd.LVectorTemp = new double[my_emd.dim];
	my_emd.RVectorTemp = new double[my_emd.dim];
}

double UB_hilbert(double*q, double*p, emd_struct& my_emd)
{
	int leftIndex = 0;
	int rightIndex = 0;
	int IndexSum = 0;
	bool nextLeave = false;

	double UB_Value = 0;

	int qMapIndex;
	int pMapIndex;

	//assignment from q and p to qTempVector and pTempVector
	for (int d = 0; d < my_emd.dim; d++)
	{
		my_emd.LVectorTemp[d] = q[d];
		my_emd.RVectorTemp[d] = p[d];
	}

	while (nextLeave == false)
	{
		if (leftIndex == my_emd.dim || rightIndex == my_emd.dim)
			break;

		if (leftIndex == my_emd.dim - 1 && rightIndex == my_emd.dim - 1)
			nextLeave = true;

		qMapIndex = my_emd.hilbert_DimVector[leftIndex];
		pMapIndex = my_emd.hilbert_DimVector[rightIndex];

		if (my_emd.LVectorTemp[qMapIndex] < my_emd.RVectorTemp[pMapIndex])
		{
			UB_Value = UB_Value + my_emd.costMatrix[qMapIndex][pMapIndex] * my_emd.LVectorTemp[qMapIndex];
			my_emd.RVectorTemp[pMapIndex] -= my_emd.LVectorTemp[qMapIndex];
			my_emd.LVectorTemp[qMapIndex] = 0;

			leftIndex++;
		}
		if (my_emd.RVectorTemp[pMapIndex] < my_emd.LVectorTemp[qMapIndex])
		{
			UB_Value = UB_Value + my_emd.costMatrix[pMapIndex][qMapIndex] * my_emd.RVectorTemp[pMapIndex];
			my_emd.LVectorTemp[qMapIndex] -= my_emd.RVectorTemp[pMapIndex];
			my_emd.RVectorTemp[pMapIndex] = 0;

			rightIndex++;
		}

		if (my_emd.LVectorTemp[qMapIndex] == my_emd.RVectorTemp[pMapIndex])
		{
			if (leftIndex == rightIndex)
			{
				leftIndex++;
				rightIndex++;
			}
			else
			{
				if (leftIndex < rightIndex)
					leftIndex++;
				else
					rightIndex++;
			}
		}

	}
	return UB_Value;
}

void UB_Hilbert_Pair_Value_Calculation(emd_struct& my_emd)
{
	for (int p = 0; p < my_emd.pairNum; p++)
		my_emd.result_pairVector.push_back(UB_hilbert(my_emd.dataSetVector[my_emd.Pair_Array[p].L_Index], my_emd.dataSetVector[my_emd.Pair_Array[p].R_Index], my_emd));
}