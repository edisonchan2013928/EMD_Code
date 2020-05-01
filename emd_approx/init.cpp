#include "init.h"

void Pair_Array_Init(fstream& pairFile, emd_struct& my_emd)
{
	my_emd.Pair_Array = new Pair[my_emd.pairNum];
	int tempIndex;

	for (int p = 0; p < my_emd.pairNum; p++)
	{
		pairFile >> tempIndex;
		my_emd.Pair_Array[p].L_Index = tempIndex;
		pairFile >> tempIndex;
		my_emd.Pair_Array[p].R_Index = tempIndex;
	}
}

void loadObject(fstream& objectSetFile, double**& objectSetVector, int objectNum, int dim)
{
	objectSetVector = new double*[objectNum];
	for (int i = 0; i < objectNum; i++)
		objectSetVector[i] = new double[dim];

	for (int i = 0; i < objectNum; i++)
		for (int d = 0; d < dim; d++)
			objectSetFile >> objectSetVector[i][d];
}

void loadData(char*dataFileName, char*pairFileName, bool isPair, emd_struct& my_emd)
{
	fstream pairFile;
	fstream dataSetFile;

	if (isPair == true)//Pair EMD
	{
		pairFile.open(pairFileName);
		if (pairFile.is_open() == false)
		{
			cout << "Cannot Open Pair File! " << endl;
			exit(0);
		}
		pairFile >> my_emd.pairNum;
		Pair_Array_Init(pairFile, my_emd);

		pairFile.close();
	}

	dataSetFile.open(dataFileName);
	if (dataSetFile.is_open() != true)
	{
		cout << "Cannot Open File!" << endl;
		exit(0);
	}

	dataSetFile >> my_emd.dataNum;
	loadObject(dataSetFile, my_emd.dataSetVector, my_emd.dataNum, my_emd.dim);

	dataSetFile.close();
}

void createCostMatrix(char*groundFileName, emd_struct& my_emd)
{
	//double*feat;
	fstream groundFile;
	//int ground_Dim;
	//double groundValue;

	//Initialization of costMatrix
	//***************************************************//
	my_emd.costMatrix = new double*[my_emd.dim];

	for (int d = 0; d < my_emd.dim; d++)
		my_emd.costMatrix[d] = new double[my_emd.dim];
	//***************************************************//

	groundFile.open(groundFileName);
	if (groundFile.is_open() == false)
	{
		cout << "Cannot Open ground file!" << endl;
		exit(0);
	}

	groundFile >> my_emd.ground_Dim;

	my_emd.groundMatrix = new double*[my_emd.dim];
	for (int d = 0; d < my_emd.dim; d++)
		my_emd.groundMatrix[d] = new double[my_emd.ground_Dim];

	for (int d = 0; d < my_emd.dim; d++)
		for (int g = 0; g < my_emd.ground_Dim; g++)
			groundFile >> my_emd.groundMatrix[d][g];

	for (int i = 0; i < my_emd.dim; i++)
	{
		for (int j = 0; j < my_emd.dim; j++)
		{
			double tmp = 0;
			for (int k = 0; k < my_emd.ground_Dim; k++)
				tmp += pow(my_emd.groundMatrix[i][k] - my_emd.groundMatrix[j][k], 2);
				//tmp += pow(feat[i*my_emd.ground_Dim + k] - feat[j*my_emd.ground_Dim + k], 2);

			my_emd.costMatrix[i][j] = sqrt(tmp);
		}
	}

	groundFile.close();
	groundFile.clear();
}

void scaling(double**dataSetVector, emd_struct& my_emd, double scal_Factor)
{
	for (int i = 0; i < my_emd.dataNum; i++)
		for (int d = 0; d < my_emd.dim; d++)
			dataSetVector[i][d] = dataSetVector[i][d] * scal_Factor;
}

void outputResultFile(char*resultFileName, emd_struct& my_emd)
{
	fstream resultFile;
	resultFile.open(resultFileName, ios::in | ios::out | ios::trunc);

	if (resultFile.is_open() == false)
	{
		cout << "Cannot open result file!" << endl;
		exit(0);
	}

	for (int p = 0; p < my_emd.pairNum; p++)
		resultFile << my_emd.result_pairVector[p] << endl;

	resultFile.close();
}