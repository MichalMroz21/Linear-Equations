#include "Matrix.h"

#include <math.h>

Matrix::Matrix(std::size_t N, int a1, int a2, int a3) //macierz pasmowa konstruktor dla rozmiaru N i diagonali
{
	A = new double* [N];
	size = N;

	for (std::size_t i{}; i < N; i++)
		A[i] = new double[N];
	
	for (std::size_t i{}; i < N; i++)
	{
		for (std::size_t j{}; j < N; j++)
		{
			if (j == i) A[i][j] = a1; //wartosci na przekatnej - a1
			else if (j == i - 1 || j == i + 1) A[i][j] = a2; //jak na rysunku a2 i a3
			else if (j == i - 2 || j == i + 2) A[i][j] = a3;
			else A[i][j] = 0; //zera
		}
	}
}

Matrix::Matrix(const Matrix& M) //make matrix from another matrix constructor.
{
	size = M.size;
	A = new double* [size];

	for (std::size_t i{}; i < size; ++i)
		A[i] = new double[size];
	

	for (std::size_t i{}; i < size; ++i)
	{
		for (std::size_t j{}; j < size; ++j)
			A[i][j] = M.A[i][j];
		
	}
}

Matrix& Matrix::operator=(const Matrix& M) //przypisanie macierzy do macierzy
{	
	for (std::size_t i{}; i < size; ++i)
	{
		for (std::size_t j{}; j < size; ++j)
			A[i][j] = M.A[i][j];
		
	}

	return *this;
}

double* Matrix::operator*(const double* v) //mnozenie macierzy razy wektor o rozmiarze N
{
	std::size_t N = size;

	double* result = new double[N];

	for (std::size_t i{}; i < N; i++)
	{
		double sum{};

		for (std::size_t j{}; j < N; j++) sum += A[i][j] * v[j];

		result[i] = sum;
	}

	return result; //wynik to wektor o dlugosci N
}