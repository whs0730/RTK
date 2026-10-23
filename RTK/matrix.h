#pragma once

#include <string>
#include <ostream>
#include <vector>

using namespace std;

typedef vector<vector<double>> Matrix;

// 创建 r 行 c 列零矩阵。
Matrix zeros(int r, int c);

// 创建 n 阶单位矩阵。
Matrix eye(int n);

// 调试用矩阵打印。
void printMatrix(const Matrix& A, const string& name = "Matrix");

// 与课程参考TXT一致：每个矩阵行占一行，列宽14，保留6位小数；不改变计算和CSV精度。
void writeMatrix(ostream& out, const Matrix& A, const string& name);

// 矩阵加法、减法、乘法。
Matrix add(const Matrix& A, const Matrix& B);
Matrix sub(const Matrix& A, const Matrix& B);
Matrix mul(const Matrix& A, const Matrix& B);

// 矩阵转置。
Matrix transpose(const Matrix& A);

// Gauss-Jordan 求逆；奇异时返回空矩阵。
Matrix inverse(const Matrix& A);

// 标量乘矩阵。
Matrix scalar_mul(double k, const Matrix& A);
