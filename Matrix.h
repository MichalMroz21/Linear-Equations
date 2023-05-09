#pragma once
#include <iostream>

class Matrix
{
public:
	double** A;
	std::size_t size;

	Matrix(std::size_t N, int a1, int a2, int a3);
	Matrix(const Matrix& M);

	Matrix& operator=(const Matrix& M);
	double* operator*(const double* v);
};