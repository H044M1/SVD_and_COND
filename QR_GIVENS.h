#ifndef	QR_GIVENS_H
#define QR_GIVENS_H

#include "matrix.h"
#include <cmath>
#include <iomanip>

void Givens_Orthogonalization(Matrix& A, Matrix& Q, Matrix& R);
void Back_Row_Substitution(const Matrix& A, const std::vector<double>& F, std::vector<double>& RES);
#endif