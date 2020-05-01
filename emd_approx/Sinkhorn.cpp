#include "Sinkhorn.h"

void round(double**F_matrix, double*q, double*p, int dim)
{
	//initialization of X_diag
	double*X_diag = new double[dim];
	double*Y_diag = new double[dim];
	double*row_F = new double[dim];
	double*col_F_dash = new double[dim];
	double*row_F_dash_dash = new double[dim];
	double*col_F_dash_dash = new double[dim];
	double*err_r = new double[dim];
	double*err_c = new double[dim];
	double**outer_matrix = new double*[dim];
	
	for (int d = 0; d < dim; d++)
		outer_matrix[d] = new double[dim];

	row_sum(F_matrix, dim, row_F);
	for (int d = 0; d < dim; d++)
	{
		if (row_F[d] < 0.00001)
			X_diag[d] = 1;
		else
			X_diag[d] = min(q[d] / row_F[d], 1.0);
	}

	for (int r = 0; r < dim; r++)
		for (int c = 0; c < dim; c++)
			F_matrix[r][c] = F_matrix[r][c] * X_diag[r];

	col_sum(F_matrix, dim, col_F_dash);
	for (int d = 0; d < dim; d++)
	{
		if (col_F_dash[d] < 0.00001)
			Y_diag[d] = 1;
		else
			Y_diag[d] = min(p[d] / col_F_dash[d], 1.0);
	}

	for (int c = 0; c < dim; c++)
		for (int r = 0; r < dim; r++)
			F_matrix[r][c] = F_matrix[r][c] * Y_diag[c];

	row_sum(F_matrix, dim, row_F_dash_dash);
	col_sum(F_matrix, dim, col_F_dash_dash);

	for (int d = 0; d < dim; d++)
	{
		err_r[d] = q[d] - row_F_dash_dash[d];
		err_c[d] = p[d] - col_F_dash_dash[d];
	}

	outer_product(err_r, err_c, dim, outer_matrix);
	matrix_divide_constant(outer_matrix, dim, L1_norm(err_r, dim));

	for (int r = 0; r < dim; r++)
		for (int c = 0; c < dim; c++)
			F_matrix[r][c] += outer_matrix[r][c];

	delete[] X_diag;
	delete[] Y_diag;
	delete[] row_F;
	delete[] col_F_dash;
	delete[] row_F_dash_dash;
	delete[] col_F_dash_dash;
	delete[] err_r;
	delete[] err_c;
	for (int d = 0; d < dim; d++)
		delete[] outer_matrix[d];
	delete[] outer_matrix;
}

double compute_dist(double*r_A, double*c_A, double*q, double*p, int dim)
{
	double dist = 0;
	for (int d = 0; d < dim; d++)
		dist += fabs(r_A[d] - q[d]) + fabs(c_A[d] - p[d]);

	return dist;
}

double rho_func(double a, double b)
{
	if (a < 0.000001)
		return 0;

	return ((b - a) + a * log(a / b));
}

int find_largest_Index(double*vec, double*vec2, int dim)
{
	int index = -1;
	double largest_value = -inf;
	double rho_value;
	for (int d = 0; d < dim; d++)
	{
		rho_value = rho_func(vec[d], vec2[d]);
		if (rho_value > largest_value)
		{
			largest_value = rho_value;
			index = d;
		}

		//debug
		//if (largest_value > 10000 || largest_value < -10000)
		//	cout << largest_value << endl;
	}
	return index;
}

void greenkhorn(double**A, double*q, double*p, int dim, double epsilon_dash)
{
	//double*x = new double[dim];
	//double*y = new double[dim];
	double x, y;
	double*r_A = new double[dim];
	double*c_A = new double[dim];
	double A_1_norm = matrix_norm_L_1(A, dim);
	int I, J;
	double dist;
	const int num_update = 1000;
	int iter = 0;

	matrix_divide_constant(A, dim, A_1_norm);
	/*for (int d = 0; d < dim; d++)
	{
		x[d] = 0;
		y[d] = 0;
	}*/
	
	row_sum(A, dim, r_A);
	col_sum(A, dim, c_A);
	dist = compute_dist(r_A, c_A, q, p, dim);
	
	//initialization
	/*for (int i = 0; i < dim; i++)
		for (int j = 0; j < dim; j++)
			B[i][j] = A[i][j];*/

	//debug
	//static int debug_iter = 0;

	while (dist > epsilon_dash && iter < num_update)
	{
		//if (iter % 50 == 0)
		//	cout << iter << ": " << dist << endl;
		//debug
		//if (debug_iter == 0)
		//	cout << dist << endl;

		I = find_largest_Index(q, r_A, dim);
		J = find_largest_Index(p, c_A, dim);
		if (rho_func(q[I], r_A[I]) > rho_func(p[J], c_A[J]))
		{
			if (r_A[I] < 0.000001)
				x = 0;
			else
				x = q[I] / r_A[I];
			dist -= fabs(q[I] - r_A[I]);
			for (int d = 0; d < dim; d++)
				A[I][d] *= x;
			r_A[I] *= x; //Update r_A[I]

			dist += fabs(q[I] - r_A[I]);
		}
		else
		{
			if (c_A[J] < 0.000001)
				y = 0;
			else
				y = p[J] / c_A[J];
			dist -= fabs(p[J] - c_A[J]);
			for (int d = 0; d < dim; d++)
				A[d][J] *= y;
			c_A[J] *= y; //Update c_A[J]

			dist += fabs(p[J] - c_A[J]);
		}
		iter++;
	}

	//debug_iter++;
	//cout << debug_iter << endl;

	delete[] r_A;
	delete[] c_A;
}

void greenkhorn_Matlab(double**A, double*q, double*p, int dim, double epsilon_dash)
{
	double x, y;
	double*r_A = new double[dim];
	double*c_A = new double[dim];
	double A_1_norm;
	int I, J;
	double dist;
	const int num_update = 1000;
	int iter = 0;

	//A_1_norm = matrix_norm_L_1(A, dim);
	//matrix_divide_constant(A, dim, A_1_norm);
	row_sum(A, dim, r_A);
	col_sum(A, dim, c_A);
	dist = compute_dist(r_A, c_A, q, p, dim);

	while (dist > epsilon_dash && iter < num_update)
	{
		I = find_largest_Index(q, r_A, dim);
		J = find_largest_Index(p, c_A, dim);
		if (rho_func(q[I], r_A[I]) > rho_func(p[J], c_A[J]))
		{
			x = q[I] / r_A[I];

			for (int d = 0; d < dim; d++)
				A[I][d] *= x;
		}
		else
		{
			y = p[J] / c_A[J];
			dist -= fabs(p[J] - c_A[J]);
			for (int d = 0; d < dim; d++)
				A[d][J] *= y;
		}

		A_1_norm = matrix_norm_L_1(A, dim);
		matrix_divide_constant(A, dim, A_1_norm);
		row_sum(A, dim, r_A);
		col_sum(A, dim, c_A);
		dist = compute_dist(r_A, c_A, q, p, dim);

		iter++;
	}

	delete[] r_A;
	delete[] c_A;
}

void store_feature(double*q, double*p, int dim)
{
	fstream file;
	file.open("debug_vec_pair.txt", ios::in | ios::out | ios::trunc);
	if (file.is_open() == false)
	{
		cout << "Cannot open file!" << endl;
		exit(0);
	}

	for (int d = 0; d < dim; d++)
		file << q[d] << " ";
	file << endl;

	for (int d = 0; d < dim; d++)
		file << p[d] << " ";
	file << endl;

	file.close();
}

void store_output(double**P, int dim, int matrix_type)
{
	fstream file;

	if (matrix_type == 0)
		file.open("debug_expMatrix.txt", ios::in | ios::out | ios::trunc);
	if (matrix_type == 1)
		file.open("debug_horn.txt", ios::in | ios::out | ios::trunc);
	if (matrix_type == 2)
		file.open("debug_round.txt", ios::in | ios::out | ios::trunc);
	if (matrix_type == 3)
		file.open("cost_matrix.txt", ios::in | ios::out | ios::trunc);

	if (file.is_open() == false)
	{
		cout << "Cannot open file!" << endl;
		exit(0);
	}

	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
			file << P[i][j] << " ";
		file << endl;
	}

	file.close();
}

double approxOT(double**costMatrix, double*q, double*p, double abs_epsilon, int dim)
{
	double eta;
	double abs_epsilon_dash;
	double**A;
	double C_infty; //specify this C_infty
	double approx_emdValue;

	A = new double*[dim];
	for (int d = 0; d < dim; d++)
		A[d] = new double[dim];

	C_infty = matrix_norm_L_infty(costMatrix, dim);
	eta = 4 * log((double)dim) / abs_epsilon;
	abs_epsilon_dash = abs_epsilon / (8 * C_infty);

	exp_Matrix(costMatrix, A, eta, dim);

	//output cost_Matrix and exp_Matrix results
	//store_output(A, dim, 0);
	//store_output(costMatrix, dim, 3);

	//Projection method
	greenkhorn(A, q, p, dim, abs_epsilon_dash);
	//greenkhorn_Matlab(A, q, p, dim, abs_epsilon_dash);

	//debug: store q,p pair and output horn result
	//store_feature(q, p, dim);
	//store_output(A, dim, 1);

	//Round method
	round(A, q, p, dim);

	//debug: output round result
	//store_output(A, dim, 2);

	//exit(0);

	approx_emdValue = matrix_inner_product(A, costMatrix, dim);

	for (int d = 0; d < dim; d++)
		delete[] A[d];
	delete[] A;

	return approx_emdValue;
}

void optimal_transport_Sinkhorn(emd_struct& my_emd)
{
	double**dataSetVector = my_emd.dataSetVector;
	for (int p = 0; p < my_emd.pairNum; p++)
		my_emd.result_pairVector.push_back(approxOT(my_emd.costMatrix, dataSetVector[my_emd.Pair_Array[p].L_Index], dataSetVector[my_emd.Pair_Array[p].R_Index], my_emd.abs_epsilon, my_emd.dim));
}