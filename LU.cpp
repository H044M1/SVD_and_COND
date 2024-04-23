#include "matrix.h"
#include <cmath>
#include <iomanip>
#include "LU.h"

using namespace std;

int findLeadingElement(const Matrix& LU, int startRow) {
    int leadingIndex = startRow;
    for (int j = startRow + 1; j < LU.rows; j++) {
        if (abs(LU.data[j][startRow]) > abs(LU.data[leadingIndex][startRow])) {
            leadingIndex = j;
        }
    }
    if (abs(LU.data[leadingIndex][startRow]) < 1e-20) {
        throw "Gauss_Method: degenerate matrix...";
    }
    return leadingIndex;
}

pair <Matrix, vector<int>> luDecomposition(const Matrix& A) {
    Matrix LU(A.rows, A.cols);
    LU = A.copy();
    vector<int> P;
    P.resize(A.rows);
    for (int i = 0; i < P.size(); i++) P[i] = i;
    for (int i = 0; i < A.rows; i++) {
        //Finding the leading element
        int idx = findLeadingElement(LU, i);



        if (idx != i) {
            swap(LU.data[idx], LU.data[i]);
            swap(P[idx], P[i]);
        }

        for (int j = i + 1; j < LU.rows; j++) {
            double help = LU.data[j][i] / LU.data[i][i];
            LU.data[j][i] = 0;
            for (int k = i + 1; k < LU.rows; k++) {
                LU.data[j][k] -= help * LU.data[i][k];
            }
        }
    }
    for (int i = 0; i < A.rows; i++) {
        for (int j = 0; j < i; j++) {
            double sum_LikUkj = 0;
            for (int k = 0; k < j; k++) {
                sum_LikUkj += LU.data[i][k] * LU.data[k][j];
            }
            LU.data[i][j] = (A.data[(int)P[i]][j] - sum_LikUkj) / LU.data[j][j];
        }
    }
    return make_pair(LU, P);
}

vector<double> directWay(const Matrix& LU, const vector<double>& F, const vector<int>& P) {
    vector<double> RES(LU.rows);
    for (int i = 0; i < LU.rows; i++) {
        RES[i] = F[(int)P[i]];
    }
    for (int i = 0; i < LU.rows; i++) {
        for (int j = 0; j < i; j++) {
            RES[i] -= LU.data[i][j] * RES[j];
        }
    }
    return RES;
}

vector<double> backWay(const Matrix& LU, const vector<double>& F, const vector<double>& P) {
    vector<double> RES(F.size());
    for (int i = 0; i < LU.rows; i++) RES[i] = F[i];

    for (int i = F.size() - 1; i >= 0; i--) {
        if (abs(LU.data[i][i]) < 1e-20) throw "������: ������� �� 0 ";

        for (int j = i + 1; j < F.size(); j++) {
            RES[i] -= LU.data[i][j] * RES[j];
        }

        RES[i] /= LU.data[i][i];
    }
    return RES;
}