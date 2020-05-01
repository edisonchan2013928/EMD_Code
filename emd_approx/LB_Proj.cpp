#include "LB_Proj.h"

typedef vector<double> cdf_Entry;

bool recordCompare(gIndexRecord record1, gIndexRecord record2)
{
	return record1.gIndexValue < record2.gIndexValue;
}

void gIndexTableVector_Creation(emd_struct& my_emd)
{
	//initalization of gIndexTable
	vector<gIndexRecord> tempVector;
	for (int gIndex = 0; gIndex < my_emd.ground_Dim; gIndex++)
		my_emd.gIndexTableVector.push_back(tempVector);

	bool abnormal;
	int beingGG;
	gIndexRecord tempRecord;

	for (int gIndex = 0; gIndex < my_emd.ground_Dim; gIndex++)
	{
		for (int d = 0; d < my_emd.dim; d++)
		{
			abnormal = false;

			for (int gg = 0; gg < (int)my_emd.gIndexTableVector[gIndex].size(); gg++)
			{
				if (fabs(my_emd.groundMatrix[d][gIndex] - my_emd.gIndexTableVector[gIndex][gg].gIndexValue) < 0.0001)//they are equal with each other
				{
					abnormal = true;
					beingGG = gg;
					break;
				}
			}

			if (abnormal == false)
			{
				tempRecord.gIndexValue = my_emd.groundMatrix[d][gIndex];
				tempRecord.dim_Set.push_back(d);
				my_emd.gIndexTableVector[gIndex].push_back(tempRecord);
				//clear
				tempRecord.dim_Set.clear();
			}
			else
				my_emd.gIndexTableVector[gIndex][beingGG].dim_Set.push_back(d);

		}

		//Sort the record according to the increasing order of gIndexValue
		sort(my_emd.gIndexTableVector[gIndex].begin(), my_emd.gIndexTableVector[gIndex].end(), recordCompare);
	}
}

void createCDF(emd_struct& my_emd)
{
	vector< cdf_Entry > temp_cdf_Entry_Vector;
	cdf_Entry temp_Entry;
	int dimIndex;
	//initialization of dataSet_CDF
	//****************************************************************//
	for (int i = 0; i < my_emd.dataNum; i++)
		my_emd.dataSet_CDF.push_back(temp_cdf_Entry_Vector);

	for (int i = 0; i < my_emd.dataNum; i++)
		for (int gIndex = 0; gIndex < my_emd.ground_Dim; gIndex++)
			my_emd.dataSet_CDF[i].push_back(temp_Entry);

	for (int i = 0; i < my_emd.dataNum; i++)
		for (int gIndex = 0; gIndex < my_emd.ground_Dim; gIndex++)
			for (int gg = 0; gg < (int)my_emd.gIndexTableVector[gIndex].size(); gg++)
				my_emd.dataSet_CDF[i][gIndex].push_back(0);
	//****************************************************************//

	//first build-up the pdf function
	//***********************************//
	for (int i = 0; i < my_emd.dataNum; i++)
	{
		for (int gIndex = 0; gIndex < my_emd.ground_Dim; gIndex++)
		{
			for (int gg = 0; gg < (int)my_emd.gIndexTableVector[gIndex].size(); gg++)
			{
				for (int d = 0; d < (int)my_emd.gIndexTableVector[gIndex][gg].dim_Set.size(); d++)
				{
					dimIndex = my_emd.gIndexTableVector[gIndex][gg].dim_Set[d];
					my_emd.dataSet_CDF[i][gIndex][gg] += my_emd.dataSetVector[i][dimIndex];
				}
			}
		}
	}
	//***********************************//

	//Make it to be the cdf
	//***********************************//
	for (int i = 0; i < my_emd.dataNum; i++)
		for (int gIndex = 0; gIndex < my_emd.ground_Dim; gIndex++)
			for (int gg = 1; gg < (int)my_emd.gIndexTableVector[gIndex].size(); gg++)
				my_emd.dataSet_CDF[i][gIndex][gg] = my_emd.dataSet_CDF[i][gIndex][gg - 1] + my_emd.dataSet_CDF[i][gIndex][gg];
	//***********************************//
}

void projectionLB(emd_struct& my_emd)
{
	double scale_Value = sqrt((double)my_emd.ground_Dim) / my_emd.ground_Dim;
	double LB_Value;
	double EMD_1DValue;

	double LValue;
	double RValue;
	double gDifference;

	for (int p = 0; p < my_emd.pairNum; p++)
	{
		LB_Value = 0;
		for (int gIndex = 0; gIndex < my_emd.ground_Dim; gIndex++)
		{
			EMD_1DValue = 0;
			for (int gg = 0; gg < (int)my_emd.gIndexTableVector[gIndex].size() - 1; gg++)
			{
				LValue = my_emd.dataSet_CDF[my_emd.Pair_Array[p].L_Index][gIndex][gg];
				RValue = my_emd.dataSet_CDF[my_emd.Pair_Array[p].R_Index][gIndex][gg];
				gDifference = my_emd.gIndexTableVector[gIndex][gg + 1].gIndexValue - my_emd.gIndexTableVector[gIndex][gg].gIndexValue;
				EMD_1DValue += fabs(RValue - LValue)*gDifference;
			}
			LB_Value += EMD_1DValue;
		}
		LB_Value *= scale_Value;
		my_emd.result_pairVector.push_back(LB_Value);
	}
}