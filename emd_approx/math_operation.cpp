#include "math_operation.h"

void row_sum(double**matrix, int dim, double*row_sumVec)
{
	//initialization
	for (int d = 0; d < dim; d++)
		row_sumVec[d] = 0;

	for (int r = 0; r < dim; r++)
		for (int c = 0; c < dim; c++)
			row_sumVec[r] += matrix[r][c];
}

void col_sum(double**matrix, int dim, double*col_sumVec)
{
	//initialization
	for (int d = 0; d < dim; d++)
		col_sumVec[d] = 0;

	for (int c = 0; c < dim; c++)
		for (int r = 0; r < dim; r++)
			col_sumVec[c] += matrix[r][c];
}

void outer_product(double*v1, double*v2, int dim, double**outerMatrix)
{
	for (int i = 0; i < dim; i++)
		for (int j = 0; j < dim; j++)
			outerMatrix[i][j] = v1[i] * v2[j];
}

void matrix_divide_constant(double**matrix, int dim, double constant)
{
	for (int r = 0; r < dim; r++)
		for (int c = 0; c < dim; c++)
			matrix[r][c] /= constant;
}

double L1_norm(double*vector, int dim)
{
	double value = 0;
	for (int d = 0; d < dim; d++)
		value += fabs(vector[d]);

	return value;
}

double matrix_norm_L_1(double**matrix, int dim)
{
	double value = 0;
	for (int i = 0; i < dim; i++)
		for (int j = 0; j < dim; j++)
			value += fabs(matrix[i][j]);

	return value;
}

double matrix_norm_L_infty(double**matrix, int dim)
{
	double max_value = -inf;
	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			if (matrix[i][j] > max_value)
				max_value = matrix[i][j];
		}
	}
	
	return max_value;
}

void exp_Matrix(double**C, double**A, double eta, int dim)
{
	for (int i = 0; i < dim; i++)
		for (int j = 0; j < dim; j++)
			A[i][j] = exp(-eta * C[i][j]);
}

double matrix_inner_product(double**P, double**C, int dim)
{
	double value = 0;
	for (int i = 0; i < dim; i++)
		for (int j = 0; j < dim; j++)
			value += P[i][j] * C[i][j];

	return value;
}
/*void matrix_multiply(double**matrixA, double**matrixB, double**C, int m, int n, int k)
{

}*/