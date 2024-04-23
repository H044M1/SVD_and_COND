#include <iostream>
#include <vector>
#include <iomanip>
#include <chrono>
#include <cmath>
#include "matrix.h"
#include "LU.h"
#include "QR_GIVENS.h"
#include "SVD.h"
using namespace std;

void fillMatrix(Matrix& mat, int N);
double get_error(const vector<double>& X_true, const vector<double>& X_found);


int main() {
	int dims[3] = { 5,10,20};
	int iterations = 1000;
	cout << setprecision(15);
	cout << fixed;
	for (auto dim : dims) {
		Matrix R(dim, dim);
		Matrix Q(dim, dim);
		Matrix A(dim, dim);
		vector<double> X_true(dim);
		for (int i = 0; i < dim; i++) {
			X_true[i] = 1;
		}
		fillMatrix(A, dim);

		vector<double> F(dim);
		F = A.multiplyVector(X_true);
		double avg_time_qr = 0;
		double avg_time_lu = 0;
		double avg_time_SVD = 0;
		double avg_err_qr = 0;
		double avg_err_lu = 0;
		double avg_err_SVD = 0;
		double cond = 0;

		cout << "dim: " << dim << endl;
		cout << "iterations: " << iterations << endl;
		auto start_qr = chrono::high_resolution_clock::now();
		vector<double> X_found_QR;
		for (int i = 0; i < iterations; i++) {
			Q.data.assign(dim, std::vector<double>(dim, 0.0));
			R.data.assign(dim, std::vector<double>(dim, 0.0));
			X_found_QR.clear();
			X_found_QR.resize(dim, 0.0);
			Givens_Orthogonalization(A, Q, R);
			X_found_QR = Q.multiplyTransposedVector(F);
			Back_Row_Substitution(R, X_found_QR, X_found_QR);
			auto end_qr = std::chrono::high_resolution_clock::now();
			avg_time_qr += chrono::duration<double>(end_qr - start_qr).count();
			avg_err_qr += get_error(X_true, X_found_QR);
		}
		for (int i = 0; i < iterations; i++) {
			auto start_lu = std::chrono::high_resolution_clock::now();
			vector<double> X_found_LU;
			Matrix LU(dim, dim);
			X_found_LU.reserve(dim);
			X_found_LU.resize(dim, 0.0);
			pair<Matrix, vector<int>> LU_P = luDecomposition(A);
			LU = LU_P.first;
			X_found_LU = directWay(LU, F, LU_P.second);
			X_found_LU = backWay(LU, X_found_LU, X_found_LU);
			auto end_lu = std::chrono::high_resolution_clock::now();
			avg_time_lu += std::chrono::duration<double>(end_lu - start_lu).count();
			avg_err_lu += get_error(X_true, X_found_LU);
			LU_P.second.clear();
			LU.data.clear();
		}
		for (int i = 0; i < iterations; i++) {
			auto start_SVD = std::chrono::high_resolution_clock::now();
			vector<double> X_found_SVD(dim);
			pair<vector<double>, double> X_found_SVD_and_cond = SVD_solver_and_cond(A, F);
			auto end_SVD = std::chrono::high_resolution_clock::now();
			avg_time_SVD += std::chrono::duration<double>(end_SVD - start_SVD).count();
			avg_err_SVD += get_error(X_true, X_found_SVD_and_cond.first);
			cond = X_found_SVD_and_cond.second;
		}
		cout << "LU:" << endl;
		cout << "avg t for LU: " << avg_time_lu / iterations << " seconds" << endl;
		cout << "error for LU: " << avg_err_lu / iterations << endl;

		cout << "QR_GIVENS:" << endl;
		cout << "avg t for QR_GIVENS: " << avg_time_qr / iterations << " seconds" << endl;
		cout << "error for  QR_GIVENS: " << avg_err_qr / iterations << endl;

		cout << "SVD :" << endl;
		cout << "avg t for SVD: " << avg_time_SVD / iterations << " seconds" << endl;
		cout << "error for SVD: " << avg_err_SVD / iterations << endl;
		cout << "Cond: " << cond << endl;
		cout << "---------------------------------------------------------------" << endl;
	}
	return 0;
}

void fillMatrix(Matrix& mat, int N) {
    if (N <= 0) {
        std::cerr << "Error: Invalid matrix size!" << std::endl;
        return;
    }

    mat.rows = N;
    mat.cols = N;
    mat.data.resize(mat.rows);

    for (int i = 0; i < mat.rows; ++i) {
        mat.data[i].resize(mat.cols, 0.0);
        for (int j = 0; j < mat.cols; ++j) {
			mat.data[i][j] = 1/(1 + 2.1 * i + 0.2 * j);
        }
    }
}

double get_error(const vector<double>& X_true, const vector<double>& X_found) {
	double norm_x_true = 0;
	double norm_x_diff = 0;
	int n = X_true.size();
	for (int i = 0; i < n; i++) {
		norm_x_true += pow(X_true[i], 2);
	}
	norm_x_true = sqrt(norm_x_true);
	for (int i = 0; i < n; i++) {
		norm_x_diff += pow(X_true[i] - X_found[i], 2);
	}
	norm_x_diff = sqrt(norm_x_diff);
	return norm_x_diff / norm_x_true;
}
