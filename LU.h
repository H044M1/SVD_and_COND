#ifndef LU_H
#define LU_H

#include "matrix.h"
#include <cmath>
#include <iomanip>

int findLeadingElement(const Matrix& LU, int startRow);
pair <Matrix, std::vector<int>> luDecomposition(const Matrix& A);
std::vector<double> directWay(const Matrix& LU, const std::vector<double>& F, const std::vector<int>& P);
std::vector<double> backWay(const Matrix& LU, const std::vector<double>& F, const std::vector<double>& P);
#endif