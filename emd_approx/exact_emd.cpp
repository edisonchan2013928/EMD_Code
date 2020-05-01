#include "exact_emd.h"

//Cost Matrix used in emd function
float**_COST;

float dist(feature_t *F1, feature_t *F2)
{
	return _COST[*F1][*F2];
}

void EMD_PairCalculation(emd_struct& my_emd)
{
	feature_t*f1 = new feature_t[my_emd.dim];
	feature_t*f2 = new feature_t[my_emd.dim];
	float*w1 = new float[my_emd.dim];
	float*w2 = new float[my_emd.dim];
	signature_t s1;
	signature_t s2;
	float e;

	flow_t*flow;
	int flowSize;
	//add code
	//if (my_emd.isCollectFlow == true)
	if (my_emd.method == 3)
		flow = new flow_t[2 * my_emd.dim];


	//Initialization and Update the cost matrix _COST
	_COST = new float*[my_emd.dim];
	for (int d1 = 0; d1 < my_emd.dim; d1++)
		_COST[d1] = new float[my_emd.dim];

	//Update the cost matrix _COST
	for (int d1 = 0; d1 < my_emd.dim; d1++)
		for (int d2 = 0; d2 < my_emd.dim; d2++)
			_COST[d1][d2] = (float)my_emd.costMatrix[d1][d2];

	for (int d = 0; d < my_emd.dim; d++)
	{
		f1[d] = d;
		f2[d] = d;
	}

	for (int p = 0; p < my_emd.pairNum; p++)
	{
		for (int d = 0; d < my_emd.dim; d++)
		{
			w1[d] = (float)my_emd.dataSetVector[my_emd.Pair_Array[p].L_Index][d];
			w2[d] = (float)my_emd.dataSetVector[my_emd.Pair_Array[p].R_Index][d];
		}
		s1.n = my_emd.dim, s1.Features = f1, s1.Weights = w1;
		s2.n = my_emd.dim, s2.Features = f2, s2.Weights = w2;
		//if (isCollectFlow == true)
		if (my_emd.method == 3)
		{
			e = emd(&s1, &s2, dist, flow, &flowSize);
			for (int f = 0; f < 2 * my_emd.dim; f++)
				if (flow[f].amount > 0)
					my_emd.flowMatrix[flow[f].from][flow[f].to] = my_emd.flowMatrix[flow[f].from][flow[f].to] + (double)flow[f].amount;
		}
		else
		{
			e = emd(&s1, &s2, dist, NULL, NULL);
			my_emd.result_pairVector.push_back(e);
		}
	}
}

void costMatrix_Reduction(double*q, double*p, int alpha, int beta, double**costMatrix, int dim)
{
	int index_q;
	int index_p;

	_COST = new float*[alpha];
	for (int i = 0; i < alpha; i++)
		_COST[i] = new float[beta];

	index_q = 0;
	for (int i = 0; i < dim; i++)
	{
		if (q[i] != 0)
		{
			index_p = 0;
			for (int ii = 0; ii < dim; ii++)
			{
				if (p[ii] != 0)
				{
					_COST[index_q][index_p] = (float)costMatrix[i][ii];
					index_p = index_p + 1;
				}
			}
			index_q = index_q + 1;
		}
	}

}

void IgnoreZero(double*q, double*p, vector<double>& q_bar, vector<double>& p_bar, double**costMatrix, int dim)
{
	q_bar.clear();
	p_bar.clear();

	int alpha = 0;
	int beta = 0;

	for (int d = 0; d < dim; d++)
	{
		if (q[d] != 0)
		{
			q_bar.push_back(q[d]);
			alpha = alpha + 1;
		}
		if (p[d] != 0)
		{
			p_bar.push_back(p[d]);
			beta = beta + 1;
		}
	}

	costMatrix_Reduction(q, p, alpha, beta, costMatrix, dim);
}

double emdComputation(vector<double>& q_bar, vector<double>& p_bar)
{
	feature_t*f1;
	feature_t*f2;
	float*w1;
	float*w2;
	signature_t s1;
	signature_t s2;
	float e;

	int q_DIM_Num = q_bar.size();
	int p_DIM_Num = p_bar.size();

	w1 = new float[q_DIM_Num];
	f1 = new feature_t[q_DIM_Num];
	w2 = new float[p_DIM_Num];
	f2 = new feature_t[p_DIM_Num];

	for (int d1 = 0; d1 < q_DIM_Num; d1++)
	{
		w1[d1] = (float)q_bar[d1];
		f1[d1] = d1;
	}

	for (int d2 = 0; d2 < p_DIM_Num; d2++)
	{
		w2[d2] = (float)p_bar[d2];
		f2[d2] = d2;
	}

	s1.n = q_DIM_Num, s1.Features = f1, s1.Weights = w1;
	s2.n = p_DIM_Num, s2.Features = f2, s2.Weights = w2;

	e = emd(&s1, &s2, dist, NULL, NULL);

	//result_pairVector.push_back(e);

	//Clear all storages
	delete[] w1;
	delete[] f1;
	delete[] w2;
	delete[] f2;
	clearCostMatrix(q_DIM_Num);
	return e;
}

void clearCostMatrix(int q_DIM_Num)
{
	for (int i = 0; i < q_DIM_Num; i++)
		delete[] _COST[i];
	delete[] _COST;
}