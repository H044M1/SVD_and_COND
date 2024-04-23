#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <iostream>
using namespace std;

class Matrix {
public:
    int rows;
    int cols;
    std::vector<vector<double>> data;
    // ����������� �� ���������
    Matrix() : rows(0), cols(0) {}

    // ����������� ����-�������
    Matrix(int numRows, int numCols) : rows(numRows), cols(numCols) {
        data.resize(rows);
        for (int i = 0; i < rows; ++i) {
            data[i].resize(cols, 0.0);
        }
    }

    // ��������� ������� �� ������ � ���������� ������ ��� ����� �������
    Matrix scalarMultiply(double scalar) const {
        Matrix result(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result.data[i][j] = data[i][j] * scalar;
            }
        }
        return result;
    }

    // ��������� ������� �� ������, ��������� ������������ � �������� �������
    void scalarMultiplyInPlace(double scalar) {
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                data[i][j] *= scalar;
            }
        }
    }

    // ��������� ������� �� ������
    std::vector<double> multiplyVector(const std::vector<double>& vec) const {
        std::vector<double> result(rows, 0.0);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result[i] += data[i][j] * vec[j];
            }
        }
        return result;
    }

    // ��������� ����������������� ������� �� ������
    std::vector<double> multiplyTransposedVector(const std::vector<double>& vec) const {
        std::vector<double> result(cols, 0.0);
        for (int i = 0; i < cols; ++i) {
            for (int j = 0; j < rows; ++j) {
                result[i] += data[j][i] * vec[j];
            }
        }
        return result;
    }

    // ��������� ������� �� �������
    Matrix multiplyMatrix(const Matrix& other) const {
        Matrix result(rows, other.cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < other.cols; ++j) {
                for (int k = 0; k < cols; ++k) {
                    result.data[i][j] += data[i][k] * other.data[k][j];
                }
            }
        }
        return result;
    }

    // �������� ������ ��� ��������� ������ ��� ����� ������� (���������� � �� �� �������)
    void addMatrix(const Matrix& other) {
        if (rows == other.rows && cols == other.cols) {
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < cols; ++j) {
                    data[i][j] += other.data[i][j];
                }
            }
        }
        else {
            cerr << "Error: Matrix dimensions do not match for addition!" << std::endl;
        }
    }

    // ���������������� �������
    Matrix transpose() const {
        Matrix result(cols, rows);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result.data[j][i] = data[i][j];
            }
        }
        return result;
    }

    // ����������� �������
    Matrix copy() const {
        Matrix result(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result.data[i][j] = data[i][j];
            }
        }
        return result;
    }

    void Column_Transposition(int i, int j) {
        for (int row = 0; row < cols; row++) {
            double temp = data[row][i];
            data[row][i] = data[row][j];
            data[row][j] = temp;
        }
    }

    void Size_Reduction(int size) {
        for (int col = size; col < cols; col++) {
            for (int row = 0; row < rows; row++) {
                data[row][col] = 1e-15;
            }

        }
    }

    // ����� ������� �� �������
    void print() const {
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                cout << data[i][j] << " ";
            }
            cout << endl;
        }
    }
};
#endif

