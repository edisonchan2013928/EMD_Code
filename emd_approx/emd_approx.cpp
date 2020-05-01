#include "emd_approx.h"

void run_emd_approx(char*dataFileName, char*pairFileName, char*groundFileName, char*resultFileName, emd_struct& my_emd)
{
	clock_t start_s;
	clock_t stop_s;
	int method = my_emd.method;

	loadData(dataFileName, pairFileName, true, my_emd);
	createCostMatrix(groundFileName, my_emd);

	//preprocessing time	
	if (method == 1)// LB_Proj
	{
		gIndexTableVector_Creation(my_emd);
		createCDF(my_emd);
	}
	if (method == 3)
	{
		//Update this parameter "my_emd.dim_Red"
		init_dimRed(my_emd);
		kMetroid_Algorithm(my_emd);
		CostMatrixReduction_Metroid(my_emd);
		dataSetReduction(my_emd);
	}
	if (method == 2 || method == 4 || method == 5)
	{
		my_emd.LVectorTemp = new double[my_emd.dim];
		my_emd.RVectorTemp = new double[my_emd.dim];

		if (method == 2) //LB_IM
			init_component(my_emd);
		
		if (method == 4) //UB_H
			initHilbert_File((char*)"", my_emd); //modify this part

		if (method == 5) //UB_G
		{
			my_emd.isSparse = true;
			list_PreProcess(my_emd.costMatrix, my_emd.sortCostList, my_emd.dim);
		}
	}

	//online runing time
	start_s = clock();
	if (method == 0)//exact EMD
		EMD_PairCalculation(my_emd);
	if (method == 1)// LB_Proj
		projectionLB(my_emd);
	if (method == 2) //LB_IM
		IM_Pair_Value_Calculation(my_emd);
	if (method == 3) //dimension reduction
		LB_dimReduction_Pair_Value_Calculation(my_emd);
	if (method == 4) //UB_H
		UB_Hilbert_Pair_Value_Calculation(my_emd);
	if (method == 5) //UB_G
		UB_Greedy_Pair_Value_Calculation(my_emd);

	stop_s = clock();

	outputResultFile(resultFileName, my_emd);

	cout << "Method " << method << ":" << my_emd.pairNum/((double)(stop_s - start_s) / CLOCKS_PER_SEC) << endl;
}
