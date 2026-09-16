#include "matrix.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

Matrix zeros(int r, int c)
{
    return Matrix(r, vector<double>(c, 0.0));
}

// 创建 n×n 单位矩阵。
Matrix eye(int n)
{
    Matrix I = zeros(n, n);
    for (int i = 0; i < n; i++) I[i][i] = 1.0;
    return I;
}

// 打印矩阵，用于调试。
void printMatrix(const Matrix& A, const string& name)
{
    cout << name << " =" << endl;
    for (int i = 0; i < (int)A.size(); i++) {
        for (int j = 0; j < (int)A[0].size(); j++) {
            cout << setw(12) << fixed << setprecision(6) << A[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

// 矩阵加法。
Matrix add(const Matrix& A, const Matrix& B)
{
    int r = A.size();
    int c = A[0].size();

    if (r != B.size() || c !=B[0].size()) {
        throw runtime_error("Matrix size mismatch in add()");
    }

    Matrix C = zeros(r, c);
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            C[i][j] = A[i][j] + B[i][j];
        }
    }
    return C;
}

// 矩阵减法。
Matrix sub(const Matrix& A, const Matrix& B)
{
    int r = (int)A.size();
    int c = (int)A[0].size();

    if (r != (int)B.size() || c != (int)B[0].size()) {
        throw runtime_error("Matrix size mismatch in sub()");
    }

    Matrix C = zeros(r, c);
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            C[i][j] = A[i][j] - B[i][j];
        }
    }
    return C;
}

// 矩阵乘法
Matrix mul(const Matrix& A, const Matrix& B)
{
    int r1 = (int)A.size();
    int c1 = (int)A[0].size();
    int r2 = (int)B.size();
    int c2 = (int)B[0].size();

    if (c1 != r2) {
        throw runtime_error("Matrix size mismatch in mul()");
    }

    Matrix C = zeros(r1, c2);
    for (int i = 0; i < r1; i++) {
        for (int j = 0; j < c2; j++) {
            for (int k = 0; k < c1; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

// 矩阵转置
Matrix transpose(const Matrix& A)
{
    int r = (int)A.size();
    int c = (int)A[0].size();

    Matrix T = zeros(c, r);
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            T[j][i] = A[i][j];
        }
    }
    return T;
}

// 矩阵求逆：高斯-约旦消元法
Matrix inverse(const Matrix& A)
{
    int n = (int)A.size();
    if (n != (int)A[0].size()) {
        throw runtime_error("Matrix must be square in inverse()");
    }

    Matrix aug = zeros(n, 2 * n);

    // 构造增广矩阵 [A | I]。
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            aug[i][j] = A[i][j];
        }
        aug[i][i + n] = 1.0;
    }

    // 高斯-约旦消元
    for (int i = 0; i < n; i++) {
        // 选主元。
        int pivot = i;
        for (int j = i + 1; j < n; j++) {
            if (fabs(aug[j][i]) > fabs(aug[pivot][i])) {
                pivot = j;
            }
        }

        if (fabs(aug[pivot][i]) < 1e-12) {
            throw runtime_error("Matrix is singular in inverse()");
        }

        // 交换行。
        if (pivot != i) {
            swap(aug[pivot], aug[i]);
        }

        // 主元归一。
        double div = aug[i][i];
        for (int j = 0; j < 2 * n; j++) {
            aug[i][j] /= div;
        }

        // 消去其他行。
        for (int k = 0; k < n; k++) {
            if (k == i) continue;
            double factor = aug[k][i];
            for (int j = 0; j < 2 * n; j++) {
                aug[k][j] -= factor * aug[i][j];
            }
        }
    }

    // 提取右半部分，得到逆矩阵
    Matrix invA = zeros(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            invA[i][j] = aug[i][j + n];
        }
    }

    return invA;
}
// 常数乘矩阵
Matrix scalar_mul(double k,const Matrix& A)
{
    int rows = A.size();
    int cols = A[0].size();

    Matrix result(rows, vector<double>(cols, 0.0));

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            result[i][j] = k * A[i][j];
        }
    }

    return result;
}


