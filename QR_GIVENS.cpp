#include "matrix.h"
#include <cmath>
#include <iomanip>
#include "QR_GIVENS.h"

using namespace std;

const double EPS = 1e-20;

void Givens_Orthogonalization(Matrix& A, Matrix& Q, Matrix& R) {
    for (int i = 0; i < A.cols; i++) {
        Q.data[i][i] = 1.0;
    }
    double help1, help2;
    double c = 0, s = 0;

    // ������ ������� � � R
    R = A.copy();

    // �������� �������� �������: ��� ������� �������
    for (int j = 0; j < R.cols - 1; j++)
    {
        // ������������� ������ � �������
        for (int i = j + 1; i < R.rows; i++)
        {
            // ���� ��������� ������� ��� ���������� �� �������, �� ��������� ������� �������
            if (abs(R.data[i][j]) > EPS)
            {
                help1 = sqrt(pow(R.data[i][j], 2) + pow(R.data[j][j], 2));
                c = R.data[j][j] / help1;
                s = R.data[i][j] / help1;

                // A_new = Gt * A
                for (int k = j; k < R.cols; k++)
                {
                    help1 = c * R.data[j][k] + s * R.data[i][k];
                    help2 = c * R.data[i][k] - s * R.data[j][k];
                    R.data[j][k] = help1;
                    R.data[i][k] = help2;
                }

                // ����������� ������ ������� Q �� ����������������� ������� �������������� Q = Q * G
                for (int k = 0; k < Q.rows; k++) {
                    help1 = c * Q.data[k][j] + s * Q.data[k][i];
                    help2 = c * Q.data[k][i] - s * Q.data[k][j];
                    Q.data[k][j] = help1;
                    Q.data[k][i] = help2;
                }
            }
        }
    }
}

void Back_Row_Substitution(const Matrix& A, const vector<double>& F, vector<double>& RES)
{
    RES = F;

    // �������� � ��������� ������, �������� �����
    for (int i = F.size() - 1; i >= 0; i--)
    {
        if (std::abs(A.data[i][i]) < EPS)
            throw std::runtime_error("Back Row Substitution: A division by zero...");

        // ��������� �� ��������
        for (int j = i + 1; j < F.size(); j++)
        {
            RES[i] -= A.data[i][j] * RES[j];
        }

        RES[i] /= A.data[i][i];
    }
}
