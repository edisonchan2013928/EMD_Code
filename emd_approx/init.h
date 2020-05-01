#pragma once
#ifndef INIT_H
#define INIT_H

#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <fstream>
#include <time.h>
#include <stdlib.h>

using namespace std;

//Define infinity
const int intINF = 2147483647;
const double inf = 2000000000000;
const double scal_Factor = 1000;

struct Pair
{
	int L_Index;
	int R_Index;
};

//Used in LB_Proj
struct gIndexRecord
{
	double gIndexValue;
	vector<int> dim_Set;
};

//Used in LB_IM
struct component
{
	int rIndex;
	double cost;
};

//UB_G
struct index_Cost
{
	int qIndex;
	int pIndex;
	double costValue;
};

struct emd_struct
{
	int queryNum;
	int dataNum;
	int pairNum;
	double**querySetVector; 
	double**dataSetVector; 
	Pair*Pair_Array;
	double**groundMatrix;
	double**costMatrix;
	int ground_Dim;
	int dim;
	vector<double> result_pairVector;

	//Different approximation methods
	int method;

	//Used in different methods
	double*LVectorTemp;
	double*RVectorTemp;

	//Used in LB_Proj (Method=1)
	vector< vector<gIndexRecord> > gIndexTableVector;
	vector< vector< vector<double> > > dataSet_CDF;

	//Used in LB_IM (Method=2)
	vector< vector<component> > SortCostVector;

	//Used in LB_Red (Method=3)
	vector< vector<int> > center_dimList; 
	double**costMatrix_Red; 
	double**Red_dataSetVector; 
	int dim_Red;

	//Used in UB_H (Method=4)
	int*hilbert_DimVector;
	char*hilbert_FileName;

	//Used in UB_G (Method=5)
	bool isSparse;
	vector<index_Cost> sortCostList;

	//Used in dimension reduction method (Method=3)
	vector<double> q_Red_bar;
	vector<double> p_Red_bar;
	double**flowMatrix;

	//Used in Sinkhorn's method
	double abs_epsilon;
	//double**G_Matrix; //output matrix
};

void Pair_Array_Init(fstream& pairFile, emd_struct& my_emd);
void loadObject(fstream& objectSetFile, double**& objectSetVector, int objectNum, int dim);
void loadData(char*dataFileName, char*pairFileName, bool isPair, emd_struct& my_emd);
void createCostMatrix(char*groundFileName, emd_struct& my_emd); 
void scaling(double**dataSetVector, emd_struct& my_emd, double scal_Factor);
void outputResultFile(char*resultFileName, emd_struct& my_emd);

#endif